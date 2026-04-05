#!/usr/bin/env python3
"""
trim_reads.py
=============
Quality-trim a parsed Sanger read FASTA using its companion quality TSV.

Trimming strategy:
  1. Always remove --trim-ends bases from both ends (primer/dye artefacts).
  2. Apply a sliding-window quality scan from both ends:
     - scan from 5' until the first base with Q >= --min-quality
     - scan from 3' until the first base with Q >= --min-quality
  3. If the remaining sequence is shorter than --min-length, the read is
     marked as FAIL and an empty FASTA record is written (with a FAIL header).
  4. IUPAC ambiguity positions within the kept region are preserved as-is but
     flagged in the output stats TSV.

Usage:
  trim_reads.py --fasta S001_forward.fasta --quality S001_forward_quality.tsv \\
                --sample-id S001 --direction forward \\
                --min-quality 20 --min-length 200 --trim-ends 20 \\
                --out-fasta S001_forward_trimmed.fasta --out-stats S001_forward_trim_stats.tsv
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Quality-trim a Sanger read FASTA.")
    p.add_argument("--fasta",       required=True, help="Input FASTA (single record)")
    p.add_argument("--quality",     required=True, help="Quality TSV from parse_abi.py")
    p.add_argument("--sample-id",   required=True, dest="sample_id")
    p.add_argument("--direction",   required=True, choices=["forward", "reverse"])
    p.add_argument("--min-quality", type=int, default=20,  dest="min_quality")
    p.add_argument("--min-length",  type=int, default=200, dest="min_length")
    p.add_argument("--trim-ends",   type=int, default=20,  dest="trim_ends")
    p.add_argument("--out-fasta",   required=True, dest="out_fasta")
    p.add_argument("--out-stats",   required=True, dest="out_stats")
    return p.parse_args()


def read_fasta(path: str) -> tuple[str, str]:
    """Return (header_without_>, sequence) for the first record in a FASTA file."""
    header = ""
    seq_parts = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if not header:
                    header = line[1:]
            else:
                seq_parts.append(line.upper())
    return header, "".join(seq_parts)


def read_quality_tsv(path: str) -> list[int]:
    """Return list of Phred quality scores (0-based indexed)."""
    quals = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                quals.append(int(row["phred_quality"]))
            except (KeyError, ValueError):
                quals.append(0)
    return quals


def find_trim_positions(quals: list[int], fixed_trim: int, min_q: int) -> tuple[int, int]:
    """
    Determine (left, right) trim positions (0-based, right is exclusive).

    1. Apply fixed_trim from each end.
    2. Advance from fixed_trim boundary until Q >= min_q.
    """
    n = len(quals)
    left = min(fixed_trim, n)
    right = max(n - fixed_trim, left)

    # Advance left boundary past low-quality bases
    while left < right and quals[left] < min_q:
        left += 1

    # Retreat right boundary past low-quality bases
    while right > left and quals[right - 1] < min_q:
        right -= 1

    return left, right


def count_ambiguous(seq: str) -> int:
    iupac = set("RYSWKMBDHVN")
    return sum(1 for b in seq if b in iupac)


def main():
    args = parse_args()

    header, seq = read_fasta(args.fasta)
    quals = read_quality_tsv(args.quality)

    # Align quals to sequence length (pad/truncate if mismatch)
    if len(quals) < len(seq):
        quals = quals + [0] * (len(seq) - len(quals))
    else:
        quals = quals[: len(seq)]

    left, right = find_trim_positions(quals, args.trim_ends, args.min_quality)
    trimmed_seq = seq[left:right]
    trimmed_len = len(trimmed_seq)
    trimmed_quals = quals[left:right]
    mean_q = sum(trimmed_quals) / trimmed_len if trimmed_len else 0.0
    n_ambig = count_ambiguous(trimmed_seq)

    qc_pass = trimmed_len >= args.min_length
    status = "PASS" if qc_pass else "FAIL"

    # Write trimmed FASTA
    with open(args.out_fasta, "w", encoding="utf-8") as fh:
        fa_header = f">{args.sample_id}_{args.direction}_trimmed status={status}"
        if qc_pass:
            fh.write(f"{fa_header}\n{trimmed_seq}\n")
        else:
            # Write empty sequence to preserve file but mark as failed
            fh.write(
                f">{args.sample_id}_{args.direction}_trimmed "
                f"status=FAIL reason=too_short "
                f"trimmed_len={trimmed_len} min_len={args.min_length}\n\n"
            )

    # Write stats TSV (one row per read)
    with open(args.out_stats, "w", encoding="utf-8") as fh:
        fh.write(
            "sample_id\tdirection\traw_length\ttrimmed_length\t"
            "trim_left\ttrim_right\tmean_quality\tambiguous_bases\tqc_status\n"
        )
        fh.write(
            f"{args.sample_id}\t{args.direction}\t{len(seq)}\t{trimmed_len}\t"
            f"{left}\t{right}\t{mean_q:.2f}\t{n_ambig}\t{status}\n"
        )

    print(
        f"[trim_reads] {args.sample_id}_{args.direction}: "
        f"raw={len(seq)} bp → trimmed={trimmed_len} bp (Q>={args.min_quality}) [{status}]",
        file=sys.stderr,
    )

    if not qc_pass:
        print(
            f"[trim_reads] WARN: {args.sample_id}_{args.direction} failed QC "
            f"(trimmed_len={trimmed_len} < min_length={args.min_length})",
            file=sys.stderr,
        )
        # Exit 0 so the pipeline can continue; BUILD_CONSENSUS will handle empty reads


if __name__ == "__main__":
    main()
