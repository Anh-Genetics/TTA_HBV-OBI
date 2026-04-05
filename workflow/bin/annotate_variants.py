#!/usr/bin/env python3
"""
annotate_variants.py
====================
Annotate variants in an HBV consensus sequence against the OBI-associated
mutation catalogue and HBV ORF annotations.

Steps
-----
1. Read the multiple sequence alignment (consensus + references) produced by MAFFT.
2. Identify the reference row (any sequence with 'REFERENCE' or 'REF' in its ID).
   Fall back to the first non-sample sequence.
3. Map each position in the alignment back to the HBsAg / S-gene coordinate system.
4. Translate the consensus S-gene (and pre-S1, pre-S2) using the expected reading frames.
5. Identify AA changes versus the reference; look them up in the OBI catalogue.
6. Output a TSV with one row per variant.

NOTE (MVP):
  - Domain coordinates (pre-S1, pre-S2, S, MHR, a-determinant) are approximate
    and based on genotype A reference NC_003977. For other genotypes the
    alignment-based remapping provides best-effort coordinates.
  - The script handles short/partial sequences gracefully but will flag them.

Usage:
  annotate_variants.py \\
      --sample-id S001 \\
      --consensus S001_consensus.fasta \\
      --alignment S001_aligned.fasta \\
      --mutation-db obi_mutation_catalogue.tsv \\
      --out-tsv S001_variants.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    from Bio import SeqIO, SeqRecord
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "[ERROR] Biopython not found. Install with: conda install -c conda-forge biopython"
    )

# ---------------------------------------------------------------------------
# HBV domain coordinates (approximate, based on NC_003977.2 genotype A)
# S-gene nt positions (1-based, within the full genome 3182 nt)
# pre-S1: 2848–3204 + 1–57  (spans origin)
# pre-S2: 3205–3204+57 / simplified below as pre-S2 start within our amplicon
# S:      155–835 (relative to preS1 start in a typical ~1245 bp amplicon)
# MHR:    S aa 99–169
# a-det:  S aa 124–147
# ---------------------------------------------------------------------------
DOMAIN_MAP = {
    # (start_nt_in_amplicon_0based, end_exclusive, label)
    # These are approximate positions within a ~1245 bp preS1–S amplicon
    "preS1":        (0,   400),
    "preS2":        (400, 550),
    "S":            (550, 1245),
    "MHR":          (850, 1057),   # approximately S codons 99–169 within amplicon
    "a_determinant":(922, 991),    # approximately S codons 124–147
}


def parse_args():
    p = argparse.ArgumentParser(description="Annotate HBV consensus variants against OBI catalogue.")
    p.add_argument("--sample-id",   required=True, dest="sample_id")
    p.add_argument("--consensus",   required=True, help="Per-sample consensus FASTA")
    p.add_argument("--alignment",   required=True, help="MAFFT alignment FASTA (consensus + refs)")
    p.add_argument("--mutation-db", required=True, dest="mutation_db")
    p.add_argument("--out-tsv",     required=True, dest="out_tsv")
    return p.parse_args()


def read_fasta_dict(path: str) -> dict[str, str]:
    """Return {seq_id: sequence_no_gaps} dict from a FASTA file."""
    records = {}
    try:
        for rec in SeqIO.parse(path, "fasta"):
            records[rec.id] = str(rec.seq).upper()
    except Exception as exc:  # noqa: BLE001
        print(f"[WARN] Could not read FASTA '{path}': {exc}", file=sys.stderr)
    return records


def read_alignment(path: str) -> dict[str, str]:
    """Return {seq_id: gapped_sequence} from a multiple alignment FASTA."""
    records = {}
    try:
        for rec in SeqIO.parse(path, "fasta"):
            records[rec.id] = str(rec.seq).upper()
    except Exception as exc:  # noqa: BLE001
        print(f"[WARN] Could not read alignment '{path}': {exc}", file=sys.stderr)
    return records


def load_mutation_catalogue(path: str) -> list[dict]:
    """
    Load the OBI mutation catalogue TSV.
    Expected columns:
      aa_position, ref_aa, alt_aa, region, mechanism, evidence_level, notes
    """
    catalogue = []
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    row["aa_position"] = int(row["aa_position"])
                except (KeyError, ValueError):
                    pass
                catalogue.append(row)
    except FileNotFoundError:
        print(f"[WARN] Mutation catalogue not found: {path}", file=sys.stderr)
    return catalogue


def find_reference_id(alignment: dict[str, str], sample_id: str) -> str | None:
    """Pick the reference sequence from the alignment (not the sample)."""
    for seq_id in alignment:
        if seq_id == sample_id:
            continue
        if any(kw in seq_id.upper() for kw in ("REF", "REFERENCE", "NC_", "AB_", "AY_", "JN_")):
            return seq_id
    # Fallback: first non-sample sequence
    for seq_id in alignment:
        if seq_id != sample_id:
            return seq_id
    return None


def degap_and_map(gapped: str) -> tuple[str, list[int]]:
    """
    Return (ungapped_seq, map_ungapped_to_alignment_pos).
    map_ungapped_to_alignment_pos[i] = alignment column index for ungapped position i.
    """
    ungapped = []
    pos_map = []
    for col_idx, base in enumerate(gapped):
        if base != "-":
            pos_map.append(col_idx)
            ungapped.append(base)
    return "".join(ungapped), pos_map


def translate_frame(nt_seq: str, frame: int = 0) -> str:
    """Translate a nucleotide string starting at *frame* offset."""
    sub = nt_seq[frame:]
    # Pad to multiple of 3
    remainder = len(sub) % 3
    if remainder:
        sub += "N" * (3 - remainder)
    try:
        return str(Seq(sub).translate(to_stop=False))
    except Exception:  # noqa: BLE001
        return ""


def annotate(
    sample_id: str,
    consensus_seq: str,
    ref_seq: str,
    catalogue: list[dict],
) -> list[dict]:
    """
    Compare consensus_seq vs ref_seq (both ungapped) and look up variants
    in the OBI catalogue.

    Returns list of variant dicts.
    """
    variants = []

    # Translate S region (approximate frame; frame=0 for the amplicon)
    # In a full preS1–S amplicon, the S ORF typically starts around nt 550
    s_start = DOMAIN_MAP["S"][0]
    cons_s = consensus_seq[s_start:] if len(consensus_seq) > s_start else ""
    ref_s  = ref_seq[s_start:]       if len(ref_seq)  > s_start else ""

    cons_aa = translate_frame(cons_s, 0) if cons_s else ""
    ref_aa  = translate_frame(ref_s,  0) if ref_s  else ""

    # Build lookup: {aa_pos: {ref_aa: catalogue_entry}}
    cat_by_pos: dict[int, list[dict]] = {}
    for entry in catalogue:
        pos = entry.get("aa_position")
        if isinstance(pos, int):
            cat_by_pos.setdefault(pos, []).append(entry)

    # Scan amino acid differences
    max_len = min(len(cons_aa), len(ref_aa))
    for aa_pos in range(max_len):
        c_aa = cons_aa[aa_pos]
        r_aa = ref_aa[aa_pos]
        if c_aa == r_aa:
            continue

        pos_1based = aa_pos + 1
        domain = _domain_label(s_start + aa_pos * 3)

        # Look up in catalogue
        cat_entries = cat_by_pos.get(pos_1based, [])
        matched = [e for e in cat_entries if e.get("alt_aa", "").upper() == c_aa]

        if matched:
            for e in matched:
                variants.append(
                    {
                        "sample_id":      sample_id,
                        "aa_position":    pos_1based,
                        "ref_aa":         r_aa,
                        "alt_aa":         c_aa,
                        "notation":       f"s{r_aa}{pos_1based}{c_aa}",
                        "domain":         domain,
                        "mechanism":      e.get("mechanism", ""),
                        "evidence_level": e.get("evidence_level", ""),
                        "notes":          e.get("notes", ""),
                        "catalogued":     "yes",
                    }
                )
        else:
            # Novel / uncatalogued variant
            variants.append(
                {
                    "sample_id":      sample_id,
                    "aa_position":    pos_1based,
                    "ref_aa":         r_aa,
                    "alt_aa":         c_aa,
                    "notation":       f"s{r_aa}{pos_1based}{c_aa}",
                    "domain":         domain,
                    "mechanism":      "unknown",
                    "evidence_level": "C",  # novel / uncatalogued
                    "notes":          "not_in_obi_catalogue",
                    "catalogued":     "no",
                }
            )

    # Check for premature stop codons in S
    stop_positions = [i + 1 for i, aa in enumerate(cons_aa) if aa == "*"]
    expected_s_len = len(ref_aa)
    premature_stops = [p for p in stop_positions if p < expected_s_len]
    for p in premature_stops:
        variants.append(
            {
                "sample_id":      sample_id,
                "aa_position":    p,
                "ref_aa":         ref_aa[p - 1] if p <= len(ref_aa) else "?",
                "alt_aa":         "*",
                "notation":       f"s{ref_aa[p-1] if p <= len(ref_aa) else '?'}{p}*",
                "domain":         _domain_label(s_start + (p - 1) * 3),
                "mechanism":      "truncation_of_HBsAg",
                "evidence_level": "A",
                "notes":          "premature_stop_codon_in_S_ORF",
                "catalogued":     "yes",
            }
        )

    return variants


def _domain_label(nt_pos_in_amplicon: int) -> str:
    # Check most-specific (smallest range) domains first
    priority_order = ["a_determinant", "MHR", "S", "preS2", "preS1"]
    for name in priority_order:
        start, end = DOMAIN_MAP[name]
        if start <= nt_pos_in_amplicon < end:
            return name
    return "outside_amplicon"


def main():
    args = parse_args()

    consensus_seqs = read_fasta_dict(args.consensus)
    alignment = read_alignment(args.alignment)
    catalogue = load_mutation_catalogue(args.mutation_db)

    # Get consensus sequence (ungapped)
    cons_seq = ""
    for sid, seq in consensus_seqs.items():
        if args.sample_id in sid or sid == args.sample_id:
            cons_seq = seq.replace("-", "")
            break
    if not cons_seq and consensus_seqs:
        cons_seq = list(consensus_seqs.values())[0].replace("-", "")

    ref_id = find_reference_id(alignment, args.sample_id)
    ref_seq = alignment.get(ref_id, "").replace("-", "") if ref_id else ""

    if not cons_seq:
        print(
            f"[annotate_variants] WARN: {args.sample_id}: empty consensus – "
            "no variants to annotate",
            file=sys.stderr,
        )
        variants = []
    elif not ref_seq:
        print(
            f"[annotate_variants] WARN: {args.sample_id}: no reference sequence "
            "in alignment – skipping annotation",
            file=sys.stderr,
        )
        variants = []
    else:
        variants = annotate(args.sample_id, cons_seq, ref_seq, catalogue)

    # Write output TSV
    fieldnames = [
        "sample_id", "aa_position", "ref_aa", "alt_aa", "notation",
        "domain", "mechanism", "evidence_level", "notes", "catalogued",
    ]
    with open(args.out_tsv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(variants)

    print(
        f"[annotate_variants] {args.sample_id}: {len(variants)} variant(s) annotated",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
