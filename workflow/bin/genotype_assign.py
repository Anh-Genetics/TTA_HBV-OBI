#!/usr/bin/env python3
"""
genotype_assign.py
==================
Parse BLAST tabular output (fmt 6) from a search against the HBV genotype
reference panel and assign genotype (A–I) to the query sample.

BLAST format 6 fields:
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

The reference FASTA headers are expected to contain the genotype tag in the
format:  >HBV_GENOTYPE_A_...  or  >GENOTYPE=A  or similar.
The script parses common naming conventions (see _extract_genotype).

Usage:
  genotype_assign.py --blast-result blast_raw.txt --sample-id S001 \\
                     --min-pct-id 90.0 --out-tsv S001_genotype.tsv
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Assign HBV genotype from BLAST results.")
    p.add_argument("--blast-result", required=True, dest="blast_result")
    p.add_argument("--sample-id",    required=True, dest="sample_id")
    p.add_argument("--min-pct-id",   type=float, default=90.0, dest="min_pct_id")
    p.add_argument("--out-tsv",      required=True, dest="out_tsv")
    return p.parse_args()


# Regex patterns to extract genotype from a BLAST subject sequence ID
_GT_PATTERNS = [
    re.compile(r"genotype[_\-=]?([A-I])", re.IGNORECASE),
    re.compile(r"\bgt[_\-]?([A-I])\b", re.IGNORECASE),
    re.compile(r"HBV[_\-]([A-I])\b"),
    re.compile(r"_([A-I])_"),      # e.g. AB_A_12345
    re.compile(r"\.([A-I])\."),    # e.g. AB.A.12345
    re.compile(r"[_\-/]([A-I])$"), # e.g. JN642165_A
]

_SUBGT_PATTERNS = [
    re.compile(r"subgenotype[_\-=]?([A-I][0-9a-z]+)", re.IGNORECASE),
    re.compile(r"\bsub[_\-]?gt[_\-]?([A-I][0-9]+)\b", re.IGNORECASE),
    re.compile(r"_([A-I][0-9]+)_"),
]


def _extract_genotype(sseqid: str) -> tuple[str, str]:
    """Return (genotype, subgenotype) or ('UNKNOWN', '') if not parseable."""
    subgt = ""
    gt = "UNKNOWN"

    for pat in _SUBGT_PATTERNS:
        m = pat.search(sseqid)
        if m:
            subgt = m.group(1).upper()
            gt = subgt[0]
            return gt, subgt

    for pat in _GT_PATTERNS:
        m = pat.search(sseqid)
        if m:
            gt = m.group(1).upper()
            if gt in "ABCDEFGHI":
                return gt, subgt

    return gt, subgt


def parse_blast_fmt6(path: str) -> list[dict]:
    """Parse BLAST format 6 output into list of dicts."""
    fields = [
        "qseqid", "sseqid", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore",
    ]
    hits = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 12:
                    continue
                row = dict(zip(fields, parts))
                try:
                    row["pident"] = float(row["pident"])
                    row["length"] = int(row["length"])
                    row["bitscore"] = float(row["bitscore"])
                    row["evalue"] = float(row["evalue"])
                except ValueError:
                    continue
                hits.append(row)
    except FileNotFoundError:
        pass  # empty results → return []
    return hits


def assign_genotype(
    hits: list[dict], min_pct_id: float
) -> tuple[str, str, float, str, str]:
    """
    Assign genotype by majority vote of top-10 hits above pident threshold.

    Returns (genotype, subgenotype, best_pct_id, best_hit_id, confidence)
    where confidence is HIGH / MEDIUM / LOW.
    """
    if not hits:
        return "UNKNOWN", "", 0.0, "", "LOW"

    qualifying = [h for h in hits if h["pident"] >= min_pct_id]
    if not qualifying:
        # Report best available hit as LOW confidence
        best = max(hits, key=lambda h: h["bitscore"])
        gt, subgt = _extract_genotype(best["sseqid"])
        return gt, subgt, best["pident"], best["sseqid"], "LOW"

    # Vote on genotype by accumulating bitscore
    gt_scores: dict[str, float] = defaultdict(float)
    subgt_scores: dict[str, float] = defaultdict(float)
    best_hit = max(qualifying, key=lambda h: h["bitscore"])

    for h in qualifying:
        gt, subgt = _extract_genotype(h["sseqid"])
        gt_scores[gt] += h["bitscore"]
        if subgt:
            subgt_scores[subgt] += h["bitscore"]

    top_gt = max(gt_scores, key=lambda k: gt_scores[k])
    top_subgt = max(subgt_scores, key=lambda k: subgt_scores[k]) if subgt_scores else ""
    if top_subgt and not top_subgt.startswith(top_gt):
        top_subgt = ""  # discard inconsistent subgenotype

    total_score = sum(gt_scores.values())
    gt_fraction = gt_scores[top_gt] / total_score if total_score else 0.0
    confidence = "HIGH" if gt_fraction >= 0.8 else ("MEDIUM" if gt_fraction >= 0.5 else "LOW")

    return top_gt, top_subgt, best_hit["pident"], best_hit["sseqid"], confidence


def main():
    args = parse_args()
    hits = parse_blast_fmt6(args.blast_result)
    genotype, subgt, pct_id, best_hit, confidence = assign_genotype(hits, args.min_pct_id)

    note = ""
    if not hits:
        note = "no_blast_hits"
    elif genotype == "UNKNOWN":
        note = f"no_hit_above_pct_id_{args.min_pct_id}"
    elif confidence == "LOW":
        note = "low_confidence_genotype_assignment"

    with open(args.out_tsv, "w", encoding="utf-8") as fh:
        fh.write(
            "sample_id\tgenotype\tsubgenotype\tpct_identity\t"
            "blast_hit\tconfidence\tnote\n"
        )
        fh.write(
            f"{args.sample_id}\t{genotype}\t{subgt}\t"
            f"{pct_id:.2f}\t{best_hit}\t{confidence}\t{note}\n"
        )

    print(
        f"[genotype_assign] {args.sample_id}: genotype={genotype} ({subgt}), "
        f"pct_id={pct_id:.1f}%, confidence={confidence}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
