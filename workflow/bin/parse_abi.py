#!/usr/bin/env python3
"""
parse_abi.py
============
Parse a single Sanger ABI (.ab1) file into:
  - A FASTA file containing the base-called sequence
  - A TSV file with per-position quality (Phred) scores and IUPAC ambiguity flags

Uses Biopython's SeqIO ABI reader, which exposes the PHRED quality scores
stored in the DATA9 / PCON1 / PCON2 channels of the .ab1 file.

Usage:
  parse_abi.py --input sample.ab1 --sample-id S001 --direction forward \\
               --out-fasta S001_forward.fasta --out-qual S001_forward_quality.tsv
"""

import argparse
import sys
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit(
        "[ERROR] Biopython is not installed.\n"
        "Install with: conda install -c conda-forge biopython\n"
        "         or: pip install biopython"
    )

# IUPAC ambiguity bases (anything other than unambiguous ACGT)
IUPAC_AMBIGUOUS = set("RYSWKMBDHVN")


def parse_args():
    p = argparse.ArgumentParser(description="Parse ABI .ab1 file to FASTA + quality TSV.")
    p.add_argument("--input",      required=True, help="Input .ab1 file")
    p.add_argument("--sample-id",  required=True, dest="sample_id")
    p.add_argument("--direction",  required=True, choices=["forward", "reverse"])
    p.add_argument("--out-fasta",  required=True, dest="out_fasta")
    p.add_argument("--out-qual",   required=True, dest="out_qual")
    return p.parse_args()


def read_abi(path: Path):
    """
    Read an ABI file with Biopython.
    Returns a SeqRecord with letter_annotations['phred_quality'].
    """
    try:
        record = SeqIO.read(str(path), "abi")
    except FileNotFoundError:
        sys.exit(f"[ERROR] ABI file not found: {path}")
    except Exception as exc:  # noqa: BLE001
        sys.exit(f"[ERROR] Failed to read ABI file '{path}': {exc}")

    # Ensure quality scores are present
    if "phred_quality" not in record.letter_annotations:
        # Try to pull from ABI-specific annotations as a fallback
        abi_data = record.annotations.get("abif_raw", {})
        pcon = abi_data.get("PCON2") or abi_data.get("PCON1")
        if pcon:
            record.letter_annotations["phred_quality"] = list(pcon[: len(record)])
        else:
            print(
                f"[WARN] No Phred quality scores found in '{path}'. "
                "Setting all qualities to 0.",
                file=sys.stderr,
            )
            record.letter_annotations["phred_quality"] = [0] * len(record)

    return record


def write_fasta(record: SeqRecord, sample_id: str, direction: str, out_path: str):
    seq_str = str(record.seq).upper()
    header = f">{sample_id}_{direction}"
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(f"{header}\n{seq_str}\n")


def write_quality_tsv(record: SeqRecord, sample_id: str, direction: str, out_path: str):
    """
    Columns:
      position (1-based), base, phred_quality, is_ambiguous
    """
    seq_str = str(record.seq).upper()
    quals = record.letter_annotations["phred_quality"]

    # Pad or trim qualities to match sequence length
    if len(quals) < len(seq_str):
        quals = list(quals) + [0] * (len(seq_str) - len(quals))
    else:
        quals = list(quals[: len(seq_str)])

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("sample_id\tdirection\tposition\tbase\tphred_quality\tis_ambiguous\n")
        for i, (base, q) in enumerate(zip(seq_str, quals), start=1):
            ambiguous = "true" if base in IUPAC_AMBIGUOUS else "false"
            fh.write(f"{sample_id}\t{direction}\t{i}\t{base}\t{q}\t{ambiguous}\n")

    # Summary statistics to stderr
    total = len(seq_str)
    mean_q = sum(quals) / total if total else 0
    n_ambig = sum(1 for b in seq_str if b in IUPAC_AMBIGUOUS)
    print(
        f"[parse_abi] {sample_id}_{direction}: {total} bp, "
        f"mean Q={mean_q:.1f}, ambiguous={n_ambig}",
        file=sys.stderr,
    )


def main():
    args = parse_args()
    ab1_path = Path(args.input)

    record = read_abi(ab1_path)
    write_fasta(record, args.sample_id, args.direction, args.out_fasta)
    write_quality_tsv(record, args.sample_id, args.direction, args.out_qual)


if __name__ == "__main__":
    main()
