#!/usr/bin/env python3
"""
trim_reads.py
=============
[EN] Quality-trim a parsed Sanger read FASTA using its companion quality TSV.
[VI] Cắt tỉa chất lượng một đoạn đọc Sanger đã phân tích dùng TSV chất lượng đi kèm.

[EN] Trimming strategy:
  1. Always remove --trim-ends bases from both ends (primer/dye artefacts).
  2. Apply a sliding-window quality scan from both ends:
     - scan from 5' until the first base with Q >= --min-quality
     - scan from 3' until the first base with Q >= --min-quality
  3. If remaining sequence < --min-length: mark as FAIL, write empty FASTA.
  4. IUPAC ambiguity positions are preserved and flagged in stats TSV.

[VI] Chiến lược cắt tỉa:
  1. Luôn bỏ --trim-ends base từ hai đầu (loại bỏ mồi và artefact thuốc nhuộm).
  2. Quét chất lượng từ hai đầu:
     - Quét từ 5' đến base đầu tiên có Q >= --min-quality
     - Quét từ 3' đến base đầu tiên có Q >= --min-quality
  3. Nếu trình tự còn lại < --min-length: đánh dấu FAIL, ghi FASTA rỗng.
  4. Các vị trí IUPAC mơ hồ được giữ nguyên và ghi lại trong TSV thống kê.

[EN] Usage:
  trim_reads.py --fasta S001_forward.fasta --quality S001_forward_quality.tsv \\
                --sample-id S001 --direction forward \\
                --min-quality 20 --min-length 200 --trim-ends 20 \\
                --out-fasta S001_forward_trimmed.fasta --out-stats S001_forward_trim_stats.tsv

[VI] Cách dùng:
  trim_reads.py --fasta M001_forward.fasta --quality M001_forward_quality.tsv \\
                --sample-id M001 --direction forward \\
                --min-quality 20 --min-length 200 --trim-ends 20 \\
                --out-fasta M001_forward_trimmed.fasta --out-stats M001_forward_trim_stats.tsv
"""

import argparse
import csv
import sys
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
            "[EN] Quality-trim a Sanger read FASTA.\n"
            "[VI] Cắt tỉa chất lượng một đoạn đọc Sanger dạng FASTA."
        )
    )
    p.add_argument("--fasta",       required=True,
                   help="[EN] Input FASTA (single record) / [VI] Tệp FASTA đầu vào")
    p.add_argument("--quality",     required=True,
                   help="[EN] Quality TSV from parse_abi.py / [VI] TSV chất lượng từ parse_abi.py")
    p.add_argument("--sample-id",   required=True, dest="sample_id",
                   help="[EN] Sample identifier / [VI] Mã mẫu")
    p.add_argument("--direction",   required=True, choices=["forward", "reverse"],
                   help="[EN] Read direction / [VI] Hướng đọc")
    p.add_argument("--min-quality", type=int, default=20,  dest="min_quality",
                   help="[EN] Minimum Phred quality / [VI] Chất lượng Phred tối thiểu (mặc định: 20)")
    p.add_argument("--min-length",  type=int, default=200, dest="min_length",
                   help="[EN] Minimum sequence length after trim / [VI] Độ dài tối thiểu sau cắt (mặc định: 200 bp)")
    p.add_argument("--trim-ends",   type=int, default=20,  dest="trim_ends",
                   help="[EN] Fixed bases to trim from each end / [VI] Số base cố định cắt mỗi đầu (mặc định: 20)")
    p.add_argument("--out-fasta",   required=True, dest="out_fasta",
                   help="[EN] Output trimmed FASTA / [VI] FASTA sau cắt tỉa")
    p.add_argument("--out-stats",   required=True, dest="out_stats",
                   help="[EN] Output trim statistics TSV / [VI] TSV thống kê cắt tỉa")
    return p.parse_args()


def read_fasta(path: str) -> tuple[str, str]:
    """
    [EN] Return (header_without_>, sequence) for the first record in a FASTA file.
    [VI] Trả về (tiêu_đề_không_có_>, trình_tự) cho bản ghi đầu tiên trong tệp FASTA.
    """
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
    """
    [EN] Return list of Phred quality scores (0-based indexed) from TSV.
    [VI] Trả về danh sách điểm chất lượng Phred (chỉ số từ 0) từ tệp TSV.
    """
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
    [EN] Determine (left, right) trim positions (0-based, right is exclusive).
         1. Apply fixed_trim from each end.
         2. Advance from fixed_trim boundary until Q >= min_q.
    [VI] Xác định vị trí (trái, phải) cắt tỉa (chỉ số 0, phải là exclusive).
         1. Áp dụng fixed_trim từ mỗi đầu.
         2. Tiến từ ranh giới fixed_trim đến khi Q >= min_q.
    """
    n = len(quals)
    # Cắt cố định từ hai đầu / Fixed trim from both ends
    left  = min(fixed_trim, n)
    right = max(n - fixed_trim, left)

    # Tiến ranh giới trái qua vùng chất lượng thấp / Advance left past low-quality
    while left < right and quals[left] < min_q:
        left += 1

    # Lùi ranh giới phải qua vùng chất lượng thấp / Retreat right past low-quality
    while right > left and quals[right - 1] < min_q:
        right -= 1

    return left, right


def count_ambiguous(seq: str) -> int:
    """
    [EN] Count IUPAC ambiguous bases in a sequence.
    [VI] Đếm số base IUPAC mơ hồ trong trình tự.
    """
    iupac = set("RYSWKMBDHVN")
    return sum(1 for b in seq if b in iupac)


def main():
    args = parse_args()

    _progress(
        f"[trim_reads] === Bắt đầu cắt tỉa / Starting trimming: "
        f"{args.sample_id} ({args.direction}) ==="
    )
    _progress(
        f"[trim_reads] Tham số / Parameters: "
        f"min_quality={args.min_quality}, min_length={args.min_length}, "
        f"trim_ends={args.trim_ends}"
    )

    # ─── Bước 1: Đọc dữ liệu vào / Step 1: Load inputs ──────────────────────
    _progress(f"[trim_reads] Đọc FASTA / Reading FASTA: {args.fasta}")
    header, seq = read_fasta(args.fasta)

    _progress(f"[trim_reads] Đọc TSV chất lượng / Reading quality TSV: {args.quality}")
    quals = read_quality_tsv(args.quality)

    # Căn chỉnh độ dài chất lượng / Align quality array length to sequence
    if len(quals) < len(seq):
        quals = quals + [0] * (len(seq) - len(quals))
    else:
        quals = quals[: len(seq)]

    _progress(f"[trim_reads] Độ dài thô / Raw length: {len(seq)} bp")

    # ─── Bước 2: Tính vị trí cắt / Step 2: Calculate trim positions ─────────
    left, right = find_trim_positions(quals, args.trim_ends, args.min_quality)
    trimmed_seq   = seq[left:right]
    trimmed_len   = len(trimmed_seq)
    trimmed_quals = quals[left:right]
    mean_q  = sum(trimmed_quals) / trimmed_len if trimmed_len else 0.0
    n_ambig = count_ambiguous(trimmed_seq)

    _progress(
        f"[trim_reads] Cắt từ vị trí / Trimmed positions: [{left}:{right}] → "
        f"{trimmed_len} bp | Q trung bình / mean Q = {mean_q:.1f} | "
        f"Base mơ hồ / ambiguous = {n_ambig}"
    )

    # ─── Bước 3: Kiểm tra QC / Step 3: QC check ─────────────────────────────
    qc_pass = trimmed_len >= args.min_length
    status  = "PASS" if qc_pass else "FAIL"
    _progress(
        f"[trim_reads] Kết quả QC / QC result: {status} "
        f"(trimmed={trimmed_len} bp, min={args.min_length} bp)"
    )

    # ─── Bước 4: Ghi FASTA cắt tỉa / Step 4: Write trimmed FASTA ───────────
    with open(args.out_fasta, "w", encoding="utf-8") as fh:
        fa_header = f">{args.sample_id}_{args.direction}_trimmed status={status}"
        if qc_pass:
            fh.write(f"{fa_header}\n{trimmed_seq}\n")
        else:
            # Ghi bản ghi rỗng để đánh dấu thất bại / Write empty record to mark failure
            fh.write(
                f">{args.sample_id}_{args.direction}_trimmed "
                f"status=FAIL reason=too_short "
                f"trimmed_len={trimmed_len} min_len={args.min_length}\n\n"
            )

    # ─── Bước 5: Ghi thống kê / Step 5: Write statistics ────────────────────
    with open(args.out_stats, "w", encoding="utf-8") as fh:
        fh.write(
            "sample_id\tdirection\traw_length\ttrimmed_length\t"
            "trim_left\ttrim_right\tmean_quality\tambiguous_bases\tqc_status\n"
        )
        fh.write(
            f"{args.sample_id}\t{args.direction}\t{len(seq)}\t{trimmed_len}\t"
            f"{left}\t{right}\t{mean_q:.2f}\t{n_ambig}\t{status}\n"
        )

    if not qc_pass:
        _progress(
            f"[trim_reads] CẢNH BÁO / WARN: {args.sample_id}_{args.direction} "
            f"không đạt QC (trimmed={trimmed_len} bp < min={args.min_length} bp). "
            f"Đoạn đọc này sẽ bị bỏ qua khi xây dựng đồng thuận / "
            f"This read will be skipped during consensus building."
        )

    _progress(
        f"[trim_reads] === Hoàn thành cắt tỉa / Trimming complete: "
        f"{args.sample_id} ({args.direction}) [{status}] ==="
    )


if __name__ == "__main__":
    main()
