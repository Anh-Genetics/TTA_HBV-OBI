#!/usr/bin/env python3
"""
build_consensus.py
==================
Build a per-sample consensus sequence from one or two trimmed Sanger reads.

Algorithm
---------
Single read:
  - Use the read as-is (after QC PASS check).

Two reads (forward + reverse):
  1. Reverse-complement the reverse read.
  2. Attempt overlap assembly by finding the longest overlapping suffix of the
     forward read that aligns to a prefix of the (RC'd) reverse read with at
     least --min-overlap bases.
  3. If a valid overlap is found, merge:
       - Take all positions unique to each read as-is.
       - For overlapping positions, resolve conflicts by --conflict policy:
           'iupac'    → use IUPAC ambiguity code for the pair of bases
           'majority' → keep the base present in the longer unique region
           'n'        → place an N
  4. If no overlap ≥ --min-overlap: fall back to the longer PASS read alone
     and flag in stats.

Outputs
-------
  --out-fasta : FASTA of the consensus sequence
  --out-stats : TSV with assembly metadata

Usage:
  build_consensus.py --sample-id S001 --reads F_trimmed.fasta R_trimmed.fasta \\
                     --min-overlap 80 --conflict iupac \\
                     --out-fasta S001_consensus.fasta --out-stats S001_consensus_stats.tsv
"""

import argparse
import sys
from pathlib import Path

# IUPAC ambiguity lookup: pair of unambiguous bases → IUPAC code
IUPAC_PAIR: dict[frozenset, str] = {
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
    frozenset({"C", "G", "T"}): "B",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "C", "G"}): "V",
}


def iupac_merge(b1: str, b2: str) -> str:
    """Return IUPAC ambiguity code for two (possibly ambiguous) bases."""
    if b1 == b2:
        return b1
    pair = frozenset({b1, b2})
    return IUPAC_PAIR.get(pair, "N")


COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn",
                            "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def parse_args():
    p = argparse.ArgumentParser(description="Build per-sample consensus from Sanger reads.")
    p.add_argument("--sample-id",   required=True, dest="sample_id")
    p.add_argument("--reads",       required=True, nargs="+", help="Trimmed FASTA files (1 or 2)")
    p.add_argument("--min-overlap", type=int, default=80, dest="min_overlap")
    p.add_argument("--conflict",    choices=["iupac", "majority", "n"], default="iupac")
    p.add_argument("--out-fasta",   required=True, dest="out_fasta")
    p.add_argument("--out-stats",   required=True, dest="out_stats")
    return p.parse_args()


def read_fasta(path: str) -> tuple[str, str]:
    """Return (header_without_>, sequence_uppercase)."""
    header = ""
    parts = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if not header:
                    header = line[1:]
            else:
                parts.append(line.upper())
    return header, "".join(parts)


def is_pass_read(header: str, seq: str) -> bool:
    """A read is valid if it has sequence and its header does NOT contain 'status=FAIL'."""
    return bool(seq) and "status=FAIL" not in header


def find_overlap(fwd: str, rev_rc: str, min_overlap: int) -> int:
    """
    Find the length of the longest suffix of *fwd* that is a prefix of *rev_rc*
    with at least min_overlap matching characters.

    Returns the overlap length, or 0 if none found.
    We use simple character comparison (not pairwise alignment) for speed.
    For production use, consider using pairwise2 / edlib.
    """
    max_overlap = min(len(fwd), len(rev_rc))
    best = 0
    for ol in range(max_overlap, min_overlap - 1, -1):
        suffix = fwd[-ol:]
        prefix = rev_rc[:ol]
        # Count mismatches allowing IUPAC ambiguity
        mismatches = sum(1 for a, b in zip(suffix, prefix) if a != b and a != "N" and b != "N")
        # Allow up to 5% mismatch in overlap region
        if mismatches / ol <= 0.05:
            best = ol
            break
    return best


def merge_overlap(fwd: str, rev_rc: str, overlap: int, conflict: str) -> str:
    """Merge two reads given a known overlap length."""
    unique_fwd = fwd[:-overlap] if overlap < len(fwd) else ""
    overlap_fwd = fwd[-overlap:] if overlap > 0 else ""
    overlap_rev = rev_rc[:overlap]
    unique_rev = rev_rc[overlap:]

    merged_overlap = []
    for bf, br in zip(overlap_fwd, overlap_rev):
        if bf == br:
            merged_overlap.append(bf)
        elif conflict == "iupac":
            merged_overlap.append(iupac_merge(bf, br))
        elif conflict == "majority":
            # Use the base from whichever end has more unique sequence
            merged_overlap.append(bf if len(unique_fwd) >= len(unique_rev) else br)
        else:  # 'n'
            merged_overlap.append("N")

    return unique_fwd + "".join(merged_overlap) + unique_rev


def main():
    args = parse_args()
    reads_info = []  # list of (header, seq, direction)

    for fasta_path in args.reads:
        hdr, seq = read_fasta(fasta_path)
        # Infer direction from filename / header
        direction = "unknown"
        fname = Path(fasta_path).stem.lower()
        if "forward" in fname or "fwd" in fname or hdr.lower().endswith("_forward_trimmed"):
            direction = "forward"
        elif "reverse" in fname or "rev" in fname or hdr.lower().endswith("_reverse_trimmed"):
            direction = "reverse"
        reads_info.append((hdr, seq, direction))

    # Filter to PASS reads only
    pass_reads = [(h, s, d) for h, s, d in reads_info if is_pass_read(h, s)]

    consensus_seq = ""
    method = "unknown"
    n_conflicts = 0
    overlap_len = 0

    if len(pass_reads) == 0:
        print(
            f"[build_consensus] WARN: {args.sample_id}: no passing reads – "
            "writing empty consensus",
            file=sys.stderr,
        )
        consensus_seq = ""
        method = "no_pass_reads"

    elif len(pass_reads) == 1:
        _, seq, direction = pass_reads[0]
        consensus_seq = seq
        method = f"single_read_{direction}"
        print(
            f"[build_consensus] {args.sample_id}: single-read consensus ({direction}), "
            f"len={len(seq)} bp",
            file=sys.stderr,
        )

    else:
        # Two (or more) reads: try to use forward + reverse
        fwd_reads = [(h, s) for h, s, d in pass_reads if d == "forward"]
        rev_reads = [(h, s) for h, s, d in pass_reads if d == "reverse"]

        if not fwd_reads:
            # No forward read, just use the first available
            _, consensus_seq, d = pass_reads[0]
            method = f"single_read_{d}"
        elif not rev_reads:
            # No reverse read
            _, consensus_seq = fwd_reads[0]
            method = "single_read_forward"
        else:
            fwd_seq = fwd_reads[0][1]
            rev_seq = rev_reads[0][1]
            rev_rc = reverse_complement(rev_seq)

            overlap_len = find_overlap(fwd_seq, rev_rc, args.min_overlap)
            if overlap_len >= args.min_overlap:
                # Count conflicting positions
                fwd_overlap = fwd_seq[-overlap_len:]
                rev_overlap = rev_rc[:overlap_len]
                n_conflicts = sum(
                    1 for a, b in zip(fwd_overlap, rev_overlap) if a != b
                )
                consensus_seq = merge_overlap(fwd_seq, rev_rc, overlap_len, args.conflict)
                method = f"paired_overlap_{overlap_len}bp"
                print(
                    f"[build_consensus] {args.sample_id}: paired consensus, "
                    f"overlap={overlap_len} bp, conflicts={n_conflicts}, "
                    f"policy={args.conflict}, len={len(consensus_seq)} bp",
                    file=sys.stderr,
                )
            else:
                # Fallback: use the longer read
                if len(fwd_seq) >= len(rev_rc):
                    consensus_seq = fwd_seq
                    method = "fallback_forward_no_overlap"
                else:
                    consensus_seq = rev_rc
                    method = "fallback_reverse_rc_no_overlap"
                print(
                    f"[build_consensus] WARN: {args.sample_id}: no overlap ≥ {args.min_overlap} bp "
                    f"found; using {method.split('_')[1]} read ({len(consensus_seq)} bp)",
                    file=sys.stderr,
                )

    # Write consensus FASTA
    with open(args.out_fasta, "w", encoding="utf-8") as fh:
        fh.write(f">{args.sample_id} method={method}\n")
        if consensus_seq:
            # Wrap at 80 chars
            for i in range(0, len(consensus_seq), 80):
                fh.write(consensus_seq[i : i + 80] + "\n")
        else:
            fh.write("\n")  # empty record

    # Write stats TSV
    with open(args.out_stats, "w", encoding="utf-8") as fh:
        fh.write(
            "sample_id\tmethod\tconsensus_length\toverlap_length\t"
            "conflict_positions\tconflict_policy\tn_pass_reads\n"
        )
        fh.write(
            f"{args.sample_id}\t{method}\t{len(consensus_seq)}\t"
            f"{overlap_len}\t{n_conflicts}\t{args.conflict}\t{len(pass_reads)}\n"
        )


if __name__ == "__main__":
    main()
