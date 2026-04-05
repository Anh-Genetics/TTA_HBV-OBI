#!/usr/bin/env python3
"""
parse_abi.py
============
[EN] Parse a single Sanger ABI (.ab1) file into:
     - A FASTA file containing the base-called sequence
     - A TSV file with per-position Phred quality scores and IUPAC ambiguity flags

[VI] Phân tích một tệp Sanger ABI (.ab1) thành:
     - Tệp FASTA chứa trình tự đọc base
     - Tệp TSV ghi điểm chất lượng Phred từng vị trí và cờ IUPAC mơ hồ

[EN] Uses Biopython's SeqIO ABI reader (DATA9 / PCON1 / PCON2 channels).
[VI] Sử dụng bộ đọc ABI của Biopython (kênh DATA9 / PCON1 / PCON2).

[EN] Usage:
  parse_abi.py --input sample.ab1 --sample-id S001 --direction forward \\
               --out-fasta S001_forward.fasta --out-qual S001_forward_quality.tsv

[VI] Cách dùng:
  parse_abi.py --input mau.ab1 --sample-id M001 --direction forward \\
               --out-fasta M001_forward.fasta --out-qual M001_forward_quality.tsv
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

# Kiểm tra Biopython / Check Biopython availability
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit(
        "[ERROR] Biopython chưa được cài đặt / Biopython is not installed.\n"
        "Cài đặt bằng / Install with:\n"
        "  conda install -c conda-forge biopython\n"
        "  or: pip install biopython"
    )

# Các base IUPAC mơ hồ (không phải A/C/G/T thuần) / IUPAC ambiguous bases
IUPAC_AMBIGUOUS = set("RYSWKMBDHVN")


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
            "[EN] Parse ABI .ab1 file to FASTA + quality TSV.\n"
            "[VI] Phân tích tệp ABI .ab1 thành FASTA và TSV chất lượng."
        )
    )
    p.add_argument("--input",      required=True, help="[EN] Input .ab1 file / [VI] Tệp .ab1 đầu vào")
    p.add_argument("--sample-id",  required=True, dest="sample_id",
                   help="[EN] Sample identifier / [VI] Mã định danh mẫu")
    p.add_argument("--direction",  required=True, choices=["forward", "reverse"],
                   help="[EN] Read direction / [VI] Hướng đọc trình tự")
    p.add_argument("--out-fasta",  required=True, dest="out_fasta",
                   help="[EN] Output FASTA file / [VI] Tệp FASTA đầu ra")
    p.add_argument("--out-qual",   required=True, dest="out_qual",
                   help="[EN] Output quality TSV / [VI] Tệp TSV chất lượng đầu ra")
    return p.parse_args()


def read_abi(path: Path) -> SeqRecord:
    """
    [EN] Read an ABI file with Biopython.
         Returns a SeqRecord with letter_annotations['phred_quality'].
    [VI] Đọc tệp ABI bằng Biopython.
         Trả về SeqRecord với letter_annotations['phred_quality'].
    """
    _progress(f"[parse_abi] Đọc tệp ABI / Reading ABI file: {path}")
    try:
        record = SeqIO.read(str(path), "abi")
    except FileNotFoundError:
        sys.exit(f"[ERROR] Không tìm thấy tệp ABI / ABI file not found: {path}")
    except Exception as exc:  # noqa: BLE001
        sys.exit(f"[ERROR] Lỗi đọc tệp ABI / Failed to read ABI file '{path}': {exc}")

    # Kiểm tra điểm chất lượng có mặt / Ensure quality scores are present
    if "phred_quality" not in record.letter_annotations:
        abi_data = record.annotations.get("abif_raw", {})
        pcon = abi_data.get("PCON2") or abi_data.get("PCON1")
        if pcon:
            record.letter_annotations["phred_quality"] = list(pcon[: len(record)])
            _progress("[parse_abi] Điểm chất lượng lấy từ kênh PCON / Quality scores loaded from PCON channel.")
        else:
            _progress(
                f"[parse_abi] CẢNH BÁO / WARN: Không tìm thấy điểm Phred trong '{path}'. "
                "Đặt tất cả bằng 0 / Setting all qualities to 0."
            )
            record.letter_annotations["phred_quality"] = [0] * len(record)

    return record


def write_fasta(record: SeqRecord, sample_id: str, direction: str, out_path: str):
    """
    [EN] Write the base-called sequence as a FASTA file.
    [VI] Ghi trình tự base thành tệp FASTA.
    """
    seq_str = str(record.seq).upper()
    header = f">{sample_id}_{direction}"
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(f"{header}\n{seq_str}\n")
    _progress(f"[parse_abi] Đã ghi FASTA / FASTA written: {out_path} ({len(seq_str)} bp)")


def write_quality_tsv(record: SeqRecord, sample_id: str, direction: str, out_path: str):
    """
    [EN] Write per-position quality scores to a TSV file.
         Columns: position (1-based), base, phred_quality, is_ambiguous
    [VI] Ghi điểm chất lượng từng vị trí vào tệp TSV.
         Cột: vị trí (bắt đầu từ 1), base, phred_quality, is_ambiguous
    """
    seq_str = str(record.seq).upper()
    quals = record.letter_annotations["phred_quality"]

    # Căn chỉnh độ dài chất lượng với trình tự / Align quality length to sequence
    if len(quals) < len(seq_str):
        quals = list(quals) + [0] * (len(seq_str) - len(quals))
    else:
        quals = list(quals[: len(seq_str)])

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("sample_id\tdirection\tposition\tbase\tphred_quality\tis_ambiguous\n")
        for i, (base, q) in enumerate(zip(seq_str, quals), start=1):
            ambiguous = "true" if base in IUPAC_AMBIGUOUS else "false"
            fh.write(f"{sample_id}\t{direction}\t{i}\t{base}\t{q}\t{ambiguous}\n")

    # Thống kê tóm tắt / Summary statistics
    total = len(seq_str)
    mean_q = sum(quals) / total if total else 0
    n_ambig = sum(1 for b in seq_str if b in IUPAC_AMBIGUOUS)

    _progress(
        f"[parse_abi] {sample_id}_{direction}: "
        f"{total} bp | Q trung bình / mean Q = {mean_q:.1f} | "
        f"Base mơ hồ / ambiguous bases = {n_ambig} | "
        f"TSV: {out_path}"
    )


def main():
    args = parse_args()
    ab1_path = Path(args.input)

    _progress(
        f"[parse_abi] === Bắt đầu phân tích ABI / Starting ABI parsing: "
        f"{args.sample_id} ({args.direction}) ==="
    )

    # Đọc tệp ABI / Read ABI file
    record = read_abi(ab1_path)

    # Ghi FASTA / Write FASTA
    write_fasta(record, args.sample_id, args.direction, args.out_fasta)

    # Ghi TSV chất lượng / Write quality TSV
    write_quality_tsv(record, args.sample_id, args.direction, args.out_qual)

    _progress(
        f"[parse_abi] === Hoàn thành phân tích ABI / ABI parsing complete: "
        f"{args.sample_id} ({args.direction}) ==="
    )


if __name__ == "__main__":
    main()
