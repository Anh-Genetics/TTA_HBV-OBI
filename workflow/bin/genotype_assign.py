#!/usr/bin/env python3
"""
genotype_assign.py
==================
[EN] Parse BLAST tabular output (fmt 6) from a search against the HBV genotype
     reference panel and assign genotype (A–I) to the query sample.
[VI] Phân tích đầu ra dạng bảng BLAST (fmt 6) từ tìm kiếm trên bộ tham chiếu
     genotype HBV và gán genotype (A–I) cho mẫu truy vấn.

[EN] BLAST format 6 fields (12 columns):
     qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
[VI] 12 cột BLAST format 6:
     ID_truy_vấn ID_tham_chiếu %_đồng_nhất độ_dài mismatch gapopen ... evalue bitscore

[EN] Reference header naming convention (recognised automatically):
     >HBV_GENOTYPE_A_NC_003977   >HBV_A_AB010291   >JN642165_A
[VI] Quy ước đặt tên tiêu đề tham chiếu (tự động nhận dạng):
     >HBV_GENOTYPE_A_NC_003977   >HBV_A_AB010291   >JN642165_A

[EN] Usage:
  genotype_assign.py --blast-result blast_raw.txt --sample-id S001 \\
                     --min-pct-id 90.0 --out-tsv S001_genotype.tsv
[VI] Cách dùng:
  genotype_assign.py --blast-result blast_raw.txt --sample-id M001 \\
                     --min-pct-id 90.0 --out-tsv M001_genotype.tsv
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def _progress(msg: str):
    """
    [EN] Print a timestamped progress message to stderr.
    [VI] In thông báo tiến trình có dấu thời gian ra stderr.
    """
    ts = datetime.now(timezone.utc).strftime("%H:%M:%S")
    print(f"  [{ts}] {msg}", file=sys.stderr)


def parse_args():
    """
    [EN] Parse command-line arguments.
    [VI] Phân tích tham số dòng lệnh.
    """
    p = argparse.ArgumentParser(
        description=(
            "[EN] Assign HBV genotype from BLAST results.\n"
            "[VI] Gán genotype HBV từ kết quả BLAST."
        )
    )
    p.add_argument("--blast-result", required=True, dest="blast_result",
                   help="[EN] BLAST fmt6 output file / [VI] Tệp đầu ra BLAST fmt6")
    p.add_argument("--sample-id",    required=True, dest="sample_id",
                   help="[EN] Sample identifier / [VI] Mã mẫu")
    p.add_argument("--min-pct-id",   type=float, default=90.0, dest="min_pct_id",
                   help="[EN] Minimum % identity to call genotype / [VI] %% đồng nhất tối thiểu (mặc định: 90.0)")
    p.add_argument("--out-tsv",      required=True, dest="out_tsv",
                   help="[EN] Output genotype TSV / [VI] Tệp TSV genotype đầu ra")
    return p.parse_args()


# ─── Biểu thức chính quy nhận dạng genotype / Regex patterns for genotype ───
_GT_PATTERNS = [
    re.compile(r"genotype[_\-=]?([A-I])", re.IGNORECASE),
    re.compile(r"\bgt[_\-]?([A-I])\b",   re.IGNORECASE),
    re.compile(r"HBV[_\-]([A-I])\b"),
    re.compile(r"_([A-I])_"),
    re.compile(r"\.([A-I])\."),
    re.compile(r"[_\-/]([A-I])$"),
]

_SUBGT_PATTERNS = [
    re.compile(r"subgenotype[_\-=]?([A-I][0-9a-z]+)", re.IGNORECASE),
    re.compile(r"\bsub[_\-]?gt[_\-]?([A-I][0-9]+)\b", re.IGNORECASE),
    re.compile(r"_([A-I][0-9]+)_"),
]


def _extract_genotype(sseqid: str) -> tuple[str, str]:
    """
    [EN] Extract (genotype, subgenotype) from a BLAST subject sequence ID.
         Returns ('UNKNOWN', '') if not parseable.
    [VI] Trích xuất (genotype, subgenotype) từ ID trình tự tham chiếu BLAST.
         Trả về ('UNKNOWN', '') nếu không thể phân tích.
    """
    subgt = ""
    gt    = "UNKNOWN"

    # Thử nhận dạng subgenotype trước / Try subgenotype first
    for pat in _SUBGT_PATTERNS:
        m = pat.search(sseqid)
        if m:
            subgt = m.group(1).upper()
            gt    = subgt[0]
            return gt, subgt

    # Sau đó thử genotype / Then try genotype
    for pat in _GT_PATTERNS:
        m = pat.search(sseqid)
        if m:
            gt = m.group(1).upper()
            if gt in "ABCDEFGHI":
                return gt, subgt

    return gt, subgt


def parse_blast_fmt6(path: str) -> list[dict]:
    """
    [EN] Parse BLAST format 6 output into list of dicts.
    [VI] Phân tích đầu ra BLAST format 6 thành danh sách dict.
    """
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
                    row["pident"]   = float(row["pident"])
                    row["length"]   = int(row["length"])
                    row["bitscore"] = float(row["bitscore"])
                    row["evalue"]   = float(row["evalue"])
                except ValueError:
                    continue
                hits.append(row)
    except FileNotFoundError:
        pass   # Không có kết quả BLAST / No BLAST results → trả về rỗng
    return hits


def assign_genotype(hits: list[dict], min_pct_id: float) -> tuple[str, str, float, str, str]:
    """
    [EN] Assign genotype by majority vote (bitscore-weighted) of top hits ≥ min_pct_id.
         Returns (genotype, subgenotype, best_pct_id, best_hit_id, confidence)
         where confidence is HIGH / MEDIUM / LOW.
    [VI] Gán genotype bằng biểu quyết đa số (theo trọng số bitscore) các hit ≥ min_pct_id.
         Trả về (genotype, subgenotype, %đồng_nhất_tốt_nhất, ID_hit_tốt_nhất, độ_tin_cậy)
         với độ_tin_cậy là HIGH / MEDIUM / LOW.
    """
    if not hits:
        return "UNKNOWN", "", 0.0, "", "LOW"

    # Lọc các hit đạt ngưỡng / Filter hits above threshold
    qualifying = [h for h in hits if h["pident"] >= min_pct_id]

    if not qualifying:
        # Báo cáo hit tốt nhất nhưng với độ tin cậy thấp / Best available hit, LOW confidence
        best = max(hits, key=lambda h: h["bitscore"])
        gt, subgt = _extract_genotype(best["sseqid"])
        return gt, subgt, best["pident"], best["sseqid"], "LOW"

    # Tích luỹ bitscore theo genotype / Accumulate bitscore per genotype
    gt_scores:    dict[str, float] = defaultdict(float)
    subgt_scores: dict[str, float] = defaultdict(float)
    best_hit = max(qualifying, key=lambda h: h["bitscore"])

    for h in qualifying:
        gt, subgt = _extract_genotype(h["sseqid"])
        gt_scores[gt] += h["bitscore"]
        if subgt:
            subgt_scores[subgt] += h["bitscore"]

    top_gt    = max(gt_scores,    key=lambda k: gt_scores[k])
    top_subgt = max(subgt_scores, key=lambda k: subgt_scores[k]) if subgt_scores else ""

    # Loại bỏ subgenotype không nhất quán / Discard inconsistent subgenotype
    if top_subgt and not top_subgt.startswith(top_gt):
        top_subgt = ""

    # Tính độ tin cậy / Calculate confidence
    total_score  = sum(gt_scores.values())
    gt_fraction  = gt_scores[top_gt] / total_score if total_score else 0.0
    confidence   = "HIGH" if gt_fraction >= 0.8 else ("MEDIUM" if gt_fraction >= 0.5 else "LOW")

    return top_gt, top_subgt, best_hit["pident"], best_hit["sseqid"], confidence


def main():
    args = parse_args()

    _progress(
        f"[genotype_assign] === Bắt đầu gán genotype / Starting genotype assignment: "
        f"{args.sample_id} | ngưỡng / threshold={args.min_pct_id}% ==="
    )

    # ─── Bước 1: Đọc kết quả BLAST / Step 1: Read BLAST results ─────────────
    _progress(f"[genotype_assign] Đọc kết quả BLAST / Reading BLAST results: {args.blast_result}")
    hits = parse_blast_fmt6(args.blast_result)
    _progress(f"[genotype_assign] Tổng số hit / Total BLAST hits: {len(hits)}")

    # ─── Bước 2: Gán genotype / Step 2: Assign genotype ─────────────────────
    _progress("[genotype_assign] Phân tích và biểu quyết genotype / Analysing and voting genotype...")
    genotype, subgt, pct_id, best_hit, confidence = assign_genotype(hits, args.min_pct_id)

    # Xây dựng ghi chú / Build note
    note = ""
    if not hits:
        note = "no_blast_hits"
    elif genotype == "UNKNOWN":
        note = f"no_hit_above_pct_id_{args.min_pct_id}"
    elif confidence == "LOW":
        note = "low_confidence_genotype_assignment"

    _progress(
        f"[genotype_assign] Kết quả / Result: genotype={genotype} ({subgt}) | "
        f"%ID={pct_id:.1f}% | độ_tin_cậy / confidence={confidence}"
    )
    if note:
        _progress(f"[genotype_assign] Ghi chú / Note: {note}")

    # ─── Bước 3: Ghi TSV / Step 3: Write TSV ────────────────────────────────
    with open(args.out_tsv, "w", encoding="utf-8") as fh:
        fh.write(
            "sample_id\tgenotype\tsubgenotype\tpct_identity\t"
            "blast_hit\tconfidence\tnote\n"
        )
        fh.write(
            f"{args.sample_id}\t{genotype}\t{subgt}\t"
            f"{pct_id:.2f}\t{best_hit}\t{confidence}\t{note}\n"
        )

    _progress(
        f"[genotype_assign] === Hoàn thành gán genotype / Genotype assignment complete: "
        f"{args.sample_id} → {genotype} [{confidence}] ==="
    )


if __name__ == "__main__":
    main()
