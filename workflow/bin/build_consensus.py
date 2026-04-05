#!/usr/bin/env python3
"""
build_consensus.py
==================
[EN] Build a per-sample consensus sequence from one (forward only) or two (F+R)
     trimmed Sanger reads.
[VI] Xây dựng trình tự đồng thuận cho từng mẫu từ một (chỉ xuôi) hoặc hai
     (xuôi+ngược) đoạn đọc Sanger đã cắt tỉa.

[EN] Algorithm:
     Single read  → use that read directly as consensus.
     Two reads    → RC the reverse read, find the longest suffix/prefix overlap
                    (≥ --min-overlap), merge with configurable conflict policy.
     No overlap   → fall back to the longer PASS read (flagged in stats).

[VI] Thuật toán:
     Một đoạn đọc → dùng nguyên đoạn đọc đó làm đồng thuận.
     Hai đoạn đọc → RC đoạn đọc ngược, tìm vùng chồng lấp suffix/prefix dài nhất
                    (≥ --min-overlap), ghép với chính sách xử lý xung đột.
     Không có vùng chồng → dùng đoạn đọc dài hơn (ghi chú trong thống kê).

[EN] Conflict policies for overlapping bases:
     'iupac'    → IUPAC ambiguity code (e.g. T+C → Y)
     'majority' → base from the longer unique region
     'n'        → place N at conflicting position

[VI] Chính sách xử lý xung đột tại vùng chồng lấp:
     'iupac'    → mã IUPAC mơ hồ (ví dụ T+C → Y)
     'majority' → base từ vùng duy nhất dài hơn
     'n'        → đặt N tại vị trí xung đột

[EN] Usage:
  build_consensus.py --sample-id S001 --reads F_trimmed.fasta R_trimmed.fasta \\
                     --min-overlap 80 --conflict iupac \\
                     --out-fasta S001_consensus.fasta --out-stats S001_consensus_stats.tsv

[VI] Cách dùng:
  build_consensus.py --sample-id M001 --reads F_cat.fasta R_cat.fasta \\
                     --min-overlap 80 --conflict iupac \\
                     --out-fasta M001_consensus.fasta --out-stats M001_consensus_stats.tsv
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

# ─── Bảng mã IUPAC cho cặp base / IUPAC ambiguity lookup for base pairs ─────
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

# Bảng bổ sung nucleotide (hỗ trợ IUPAC) / Complement table (supports IUPAC)
COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn"
)


def _progress(msg: str):
    """
    [EN] Print a timestamped progress message to stderr.
    [VI] In thông báo tiến trình có dấu thời gian ra stderr.
    """
    ts = datetime.now(timezone.utc).strftime("%H:%M:%S")
    print(f"  [{ts}] {msg}", file=sys.stderr)


def iupac_merge(b1: str, b2: str) -> str:
    """
    [EN] Return IUPAC ambiguity code for two (possibly ambiguous) bases.
    [VI] Trả về mã IUPAC mơ hồ cho hai base (có thể đã mơ hồ).
    """
    if b1 == b2:
        return b1
    pair = frozenset({b1, b2})
    return IUPAC_PAIR.get(pair, "N")


def reverse_complement(seq: str) -> str:
    """
    [EN] Return the reverse complement of a nucleotide sequence (IUPAC-aware).
    [VI] Trả về chuỗi bổ sung đảo ngược của trình tự nucleotide (hỗ trợ IUPAC).
    """
    return seq.translate(COMPLEMENT)[::-1]


def parse_args():
    """
    [EN] Parse command-line arguments.
    [VI] Phân tích tham số dòng lệnh.
    """
    p = argparse.ArgumentParser(
        description=(
            "[EN] Build per-sample consensus from Sanger reads.\n"
            "[VI] Xây dựng trình tự đồng thuận từ các đoạn đọc Sanger."
        )
    )
    p.add_argument("--sample-id",   required=True, dest="sample_id",
                   help="[EN] Sample identifier / [VI] Mã mẫu")
    p.add_argument("--reads",       required=True, nargs="+",
                   help="[EN] Trimmed FASTA files (1 or 2) / [VI] Tệp FASTA đã cắt (1 hoặc 2)")
    p.add_argument("--min-overlap", type=int, default=80, dest="min_overlap",
                   help="[EN] Min overlap to attempt merge (bp) / [VI] Vùng chồng tối thiểu để ghép (bp)")
    p.add_argument("--conflict",    choices=["iupac", "majority", "n"], default="iupac",
                   help="[EN] Conflict resolution policy / [VI] Chính sách giải quyết xung đột")
    p.add_argument("--out-fasta",   required=True, dest="out_fasta",
                   help="[EN] Output consensus FASTA / [VI] FASTA đồng thuận đầu ra")
    p.add_argument("--out-stats",   required=True, dest="out_stats",
                   help="[EN] Output consensus statistics TSV / [VI] TSV thống kê đồng thuận đầu ra")
    return p.parse_args()


def read_fasta(path: str) -> tuple[str, str]:
    """
    [EN] Return (header_without_>, sequence_uppercase).
    [VI] Trả về (tiêu_đề_không_có_>, trình_tự_in_hoa).
    """
    header, parts = "", []
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
    """
    [EN] A read passes QC if it has sequence and header does NOT contain 'status=FAIL'.
    [VI] Đoạn đọc vượt QC nếu có trình tự và tiêu đề KHÔNG chứa 'status=FAIL'.
    """
    return bool(seq) and "status=FAIL" not in header


def find_overlap(fwd: str, rev_rc: str, min_overlap: int) -> int:
    """
    [EN] Find the longest suffix of *fwd* that is also a prefix of *rev_rc*
         with ≤5% mismatches and length ≥ min_overlap.
         Returns the overlap length, or 0 if none found.
    [VI] Tìm suffix dài nhất của *fwd* cũng là prefix của *rev_rc*
         với ≤5% mismatch và độ dài ≥ min_overlap.
         Trả về độ dài vùng chồng, hoặc 0 nếu không tìm thấy.
    """
    max_overlap = min(len(fwd), len(rev_rc))
    for ol in range(max_overlap, min_overlap - 1, -1):
        suffix = fwd[-ol:]
        prefix = rev_rc[:ol]
        # Tính tỉ lệ mismatch / Calculate mismatch rate
        mismatches = sum(
            1 for a, b in zip(suffix, prefix)
            if a != b and a != "N" and b != "N"
        )
        if mismatches / ol <= 0.05:
            return ol
    return 0


def merge_overlap(fwd: str, rev_rc: str, overlap: int, conflict: str) -> str:
    """
    [EN] Merge two reads given a known overlap length using the specified conflict policy.
    [VI] Ghép hai đoạn đọc với vùng chồng lấp đã biết, dùng chính sách xử lý xung đột.
    """
    # Phần duy nhất của mỗi đoạn đọc / Unique portions of each read
    unique_fwd = fwd[:-overlap] if overlap < len(fwd) else ""
    overlap_fwd = fwd[-overlap:] if overlap > 0 else ""
    overlap_rev = rev_rc[:overlap]
    unique_rev  = rev_rc[overlap:]

    # Giải quyết xung đột tại vùng chồng / Resolve conflicts in overlap
    merged_overlap = []
    for bf, br in zip(overlap_fwd, overlap_rev):
        if bf == br:
            merged_overlap.append(bf)
        elif conflict == "iupac":
            merged_overlap.append(iupac_merge(bf, br))
        elif conflict == "majority":
            # Dùng base từ vùng duy nhất dài hơn / Use base from longer unique region
            merged_overlap.append(bf if len(unique_fwd) >= len(unique_rev) else br)
        else:   # 'n'
            merged_overlap.append("N")

    return unique_fwd + "".join(merged_overlap) + unique_rev


def main():
    args = parse_args()

    _progress(
        f"[build_consensus] === Bắt đầu xây dựng đồng thuận / Starting consensus: "
        f"{args.sample_id} | chính sách / policy={args.conflict} | "
        f"min_overlap={args.min_overlap} bp ==="
    )

    # ─── Bước 1: Đọc tất cả tệp FASTA / Step 1: Read all FASTA files ────────
    reads_info = []
    for fasta_path in args.reads:
        _progress(f"[build_consensus] Đọc / Reading: {fasta_path}")
        hdr, seq = read_fasta(fasta_path)
        # Suy ra hướng từ tên tệp / Infer direction from filename
        direction = "unknown"
        fname = Path(fasta_path).stem.lower()
        if "forward" in fname or "fwd" in fname or hdr.lower().endswith("_forward_trimmed"):
            direction = "forward"
        elif "reverse" in fname or "rev" in fname or hdr.lower().endswith("_reverse_trimmed"):
            direction = "reverse"
        reads_info.append((hdr, seq, direction))
        _progress(
            f"[build_consensus] → {direction}: {len(seq)} bp | "
            f"PASS={is_pass_read(hdr, seq)}"
        )

    # ─── Bước 2: Lọc đoạn đọc vượt QC / Step 2: Filter passing reads ────────
    pass_reads = [(h, s, d) for h, s, d in reads_info if is_pass_read(h, s)]
    _progress(
        f"[build_consensus] {len(pass_reads)}/{len(reads_info)} đoạn đọc vượt QC / "
        f"read(s) passed QC"
    )

    consensus_seq = ""
    method        = "unknown"
    n_conflicts   = 0
    overlap_len   = 0

    # ─── Bước 3: Xây dựng đồng thuận / Step 3: Build consensus ──────────────
    if len(pass_reads) == 0:
        # Không có đoạn đọc nào vượt QC / No passing reads
        _progress(
            f"[build_consensus] CẢNH BÁO / WARN: {args.sample_id}: "
            "không có đoạn đọc vượt QC / no passing reads – "
            "ghi đồng thuận rỗng / writing empty consensus"
        )
        method = "no_pass_reads"

    elif len(pass_reads) == 1:
        # Một đoạn đọc / Single read
        _, seq, direction = pass_reads[0]
        consensus_seq = seq
        method = f"single_read_{direction}"
        _progress(
            f"[build_consensus] Đồng thuận từ một đoạn đọc / Single-read consensus "
            f"({direction}): {len(seq)} bp"
        )

    else:
        # Hai đoạn đọc: thử ghép / Two reads: attempt overlap assembly
        fwd_reads = [(h, s) for h, s, d in pass_reads if d == "forward"]
        rev_reads = [(h, s) for h, s, d in pass_reads if d == "reverse"]

        if not fwd_reads:
            _, consensus_seq, d = pass_reads[0]
            method = f"single_read_{d}"
            _progress(f"[build_consensus] Không có đoạn xuôi / No forward read – dùng / using {d}")
        elif not rev_reads:
            _, consensus_seq = fwd_reads[0]
            method = "single_read_forward"
            _progress("[build_consensus] Không có đoạn ngược / No reverse read – dùng / using forward")
        else:
            fwd_seq = fwd_reads[0][1]
            rev_seq = rev_reads[0][1]
            _progress(
                f"[build_consensus] Tính bổ sung đảo ngược / RC-ing reverse read: "
                f"{len(rev_seq)} bp → RC"
            )
            rev_rc = reverse_complement(rev_seq)

            _progress(
                f"[build_consensus] Tìm vùng chồng lấp / Finding overlap "
                f"(min={args.min_overlap} bp)..."
            )
            overlap_len = find_overlap(fwd_seq, rev_rc, args.min_overlap)

            if overlap_len >= args.min_overlap:
                # Đếm vị trí xung đột / Count conflicting positions
                n_conflicts = sum(
                    1 for a, b in zip(fwd_seq[-overlap_len:], rev_rc[:overlap_len])
                    if a != b
                )
                _progress(
                    f"[build_consensus] Vùng chồng lấp / Overlap found: "
                    f"{overlap_len} bp | Xung đột / Conflicts: {n_conflicts}"
                )
                consensus_seq = merge_overlap(fwd_seq, rev_rc, overlap_len, args.conflict)
                method = f"paired_overlap_{overlap_len}bp"
                _progress(
                    f"[build_consensus] Đồng thuận cặp đọc / Paired consensus: "
                    f"{len(consensus_seq)} bp (chính sách / policy={args.conflict})"
                )
            else:
                # Không đủ chồng lấp → dùng đoạn đọc dài hơn / Insufficient overlap → use longer read
                if len(fwd_seq) >= len(rev_rc):
                    consensus_seq = fwd_seq
                    method = "fallback_forward_no_overlap"
                else:
                    consensus_seq = rev_rc
                    method = "fallback_reverse_rc_no_overlap"
                _progress(
                    f"[build_consensus] CẢNH BÁO / WARN: Không tìm thấy vùng chồng ≥ {args.min_overlap} bp. "
                    f"Dùng / Using {method.split('_')[1]} read ({len(consensus_seq)} bp)."
                )

    # ─── Bước 4: Ghi FASTA đồng thuận / Step 4: Write consensus FASTA ───────
    with open(args.out_fasta, "w", encoding="utf-8") as fh:
        fh.write(f">{args.sample_id} method={method}\n")
        if consensus_seq:
            # Ngắt dòng mỗi 80 ký tự / Wrap at 80 characters
            for i in range(0, len(consensus_seq), 80):
                fh.write(consensus_seq[i : i + 80] + "\n")
        else:
            fh.write("\n")  # Bản ghi rỗng / Empty record

    # ─── Bước 5: Ghi thống kê / Step 5: Write statistics ────────────────────
    with open(args.out_stats, "w", encoding="utf-8") as fh:
        fh.write(
            "sample_id\tmethod\tconsensus_length\toverlap_length\t"
            "conflict_positions\tconflict_policy\tn_pass_reads\n"
        )
        fh.write(
            f"{args.sample_id}\t{method}\t{len(consensus_seq)}\t"
            f"{overlap_len}\t{n_conflicts}\t{args.conflict}\t{len(pass_reads)}\n"
        )

    _progress(
        f"[build_consensus] === Hoàn thành xây dựng đồng thuận / Consensus complete: "
        f"{args.sample_id} | phương pháp / method={method} | "
        f"{len(consensus_seq)} bp ==="
    )


if __name__ == "__main__":
    main()
