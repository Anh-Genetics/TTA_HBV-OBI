#!/usr/bin/env python3
"""
validate_samplesheet.py
=======================
[EN] Validate the HBV-OBI pipeline sample sheet CSV.
[VI] Kiểm tra hợp lệ tệp mô tả mẫu (sample sheet) CSV cho pipeline HBV-OBI.

[EN] Expected CSV columns (case-insensitive, order-independent):
  sample_id     – unique sample identifier / mã định danh mẫu duy nhất
  ab1_forward   – path to forward (.ab1) read / đường dẫn đến tệp đọc xuôi
  ab1_reverse   – path to reverse (.ab1) read [OPTIONAL] / tệp đọc ngược [tuỳ chọn]

[EN] Outputs:
  --output : validated TSV with canonical column names and resolved absolute paths
  --report : human-readable validation report (text)

[VI] Đầu ra:
  --output : TSV đã kiểm tra với tên cột chuẩn và đường dẫn tuyệt đối
  --report : báo cáo kiểm tra dạng văn bản dễ đọc

[EN] Exit codes: 0 = passed, 1 = one or more fatal errors
[VI] Mã thoát: 0 = hợp lệ, 1 = có lỗi nghiêm trọng
"""

import argparse
import csv
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

# ─── Cột bắt buộc và tuỳ chọn / Required and optional columns ───────────────
REQUIRED_COLUMNS = {"sample_id", "ab1_forward"}
OPTIONAL_COLUMNS = {"ab1_reverse"}

# Định dạng tệp ABI được chấp nhận / Accepted ABI file extensions
VALID_ABI_EXTENSIONS = {".ab1", ".abi"}


def parse_args():
    """
    [EN] Parse command-line arguments.
    [VI] Phân tích các tham số dòng lệnh.
    """
    p = argparse.ArgumentParser(
        description=(
            "[EN] Validate the HBV-OBI pipeline sample sheet.\n"
            "[VI] Kiểm tra hợp lệ tệp mô tả mẫu cho pipeline HBV-OBI."
        )
    )
    p.add_argument("--input",  required=True, help="[EN] Input CSV sample sheet / [VI] Tệp CSV đầu vào")
    p.add_argument("--output", required=True, help="[EN] Output validated TSV / [VI] Tệp TSV đầu ra đã kiểm tra")
    p.add_argument("--report", required=True, help="[EN] Output validation report / [VI] Tệp báo cáo kiểm tra")
    return p.parse_args()


def normalise_header(header: list[str]) -> dict[str, str]:
    """
    [EN] Return a mapping of lower-cased column names to original names.
    [VI] Trả về ánh xạ tên cột viết thường sang tên gốc.
    """
    return {col.strip().lower(): col.strip() for col in header}


def validate_file_path(path_str: str, field: str, sample_id: str) -> tuple[str, list[str]]:
    """
    [EN] Validate that a file path exists and has a recognised ABI extension.
         Returns (resolved_absolute_path, [error_messages]).
    [VI] Kiểm tra đường dẫn tệp tồn tại và có phần mở rộng ABI hợp lệ.
         Trả về (đường_dẫn_tuyệt_đối, [danh_sách_lỗi]).
    """
    errors = []

    # Trường rỗng → tệp tuỳ chọn vắng mặt / Empty field → optional file absent
    if not path_str or path_str.strip() == "":
        return "NONE", errors

    p = Path(path_str.strip())
    # Chuyển đổi đường dẫn tương đối thành tuyệt đối / Convert relative to absolute
    if not p.is_absolute():
        p = Path.cwd() / p

    # Kiểm tra tệp tồn tại / Check file existence
    if not p.exists():
        errors.append(
            f"  [ERROR] {field} for sample '{sample_id}': "
            f"file not found / tệp không tồn tại: {p}"
        )
    # Kiểm tra phần mở rộng / Check file extension
    elif p.suffix.lower() not in VALID_ABI_EXTENSIONS:
        errors.append(
            f"  [WARN]  {field} for sample '{sample_id}': unexpected extension "
            f"'{p.suffix}' (expected .ab1 or .abi) – continuing anyway / "
            f"phần mở rộng không mong đợi – vẫn tiếp tục"
        )

    return str(p), errors


def _progress(msg: str):
    """
    [EN] Print a progress message to stderr with a timestamp.
    [VI] In thông báo tiến trình ra stderr kèm dấu thời gian.
    """
    ts = datetime.now(timezone.utc).strftime("%H:%M:%S")
    print(f"  [{ts}] {msg}", file=sys.stderr)


def main():
    args = parse_args()
    errors_total = []
    warnings_total = []
    validated_rows = []

    # ─── Bước 1: Mở và đọc CSV / Step 1: Open and read CSV ─────────────────
    _progress("[validate_samplesheet] Bắt đầu kiểm tra tệp mô tả mẫu / Starting sample sheet validation...")
    _progress(f"[validate_samplesheet] Đọc tệp đầu vào / Reading input: {args.input}")

    try:
        with open(args.input, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                sys.exit(
                    "[ERROR] Sample sheet appears to be empty / "
                    "Tệp mô tả mẫu có vẻ rỗng."
                )

            norm_map = normalise_header(list(reader.fieldnames))
            lower_cols = set(norm_map.keys())

            # ── Kiểm tra cột bắt buộc / Check required columns ──────────────
            missing = REQUIRED_COLUMNS - lower_cols
            if missing:
                sys.exit(
                    f"[ERROR] Missing required column(s) / Thiếu cột bắt buộc: "
                    f"{', '.join(sorted(missing))}.\n"
                    f"Found / Tìm thấy: {', '.join(sorted(lower_cols))}"
                )

            has_reverse = "ab1_reverse" in lower_cols
            seen_ids: set[str] = set()

            # ── Kiểm tra từng hàng mẫu / Validate each sample row ───────────
            for lineno, row in enumerate(reader, start=2):
                norm_row = {k.strip().lower(): v for k, v in row.items()}

                sample_id = norm_row.get("sample_id", "").strip()
                if not sample_id:
                    errors_total.append(
                        f"  [ERROR] Line {lineno}: empty sample_id / mã mẫu rỗng – skipping row"
                    )
                    continue

                # Kiểm tra trùng lặp mã mẫu / Check duplicate sample_id
                if sample_id in seen_ids:
                    errors_total.append(
                        f"  [ERROR] Duplicate sample_id / Trùng mã mẫu '{sample_id}' at line {lineno}"
                    )
                seen_ids.add(sample_id)

                # Kiểm tra tệp đọc xuôi / Validate forward read file
                fwd_raw = norm_row.get("ab1_forward", "")
                fwd_path, fwd_errs = validate_file_path(fwd_raw, "ab1_forward", sample_id)
                if fwd_path == "NONE":
                    errors_total.append(
                        f"  [ERROR] Line {lineno}: sample '{sample_id}' has no ab1_forward value / "
                        f"không có tệp đọc xuôi"
                    )
                else:
                    for e in fwd_errs:
                        if "[ERROR]" in e:
                            errors_total.append(e)
                        else:
                            warnings_total.append(e)

                # Kiểm tra tệp đọc ngược (tuỳ chọn) / Validate reverse read (optional)
                rev_raw = norm_row.get("ab1_reverse", "") if has_reverse else ""
                rev_path, rev_errs = validate_file_path(rev_raw, "ab1_reverse", sample_id)
                for e in rev_errs:
                    if "[ERROR]" in e:
                        errors_total.append(e)
                    else:
                        warnings_total.append(e)

                validated_rows.append(
                    {
                        "sample_id":    sample_id,
                        "ab1_forward":  fwd_path,
                        "ab1_reverse":  rev_path,
                    }
                )
                _progress(
                    f"[validate_samplesheet] Mẫu / Sample '{sample_id}': "
                    f"F={fwd_path!r}, R={rev_path!r}"
                )

    except FileNotFoundError:
        sys.exit(f"[ERROR] Sample sheet not found / Không tìm thấy tệp: {args.input}")
    except csv.Error as exc:
        sys.exit(f"[ERROR] Failed to parse CSV / Lỗi phân tích CSV: {exc}")

    # ─── Bước 2: Ghi báo cáo kiểm tra / Step 2: Write validation report ────
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    status = "PASSED" if not errors_total else "FAILED"

    report_lines = [
        "HBV-OBI Sample Sheet Validation Report / Báo cáo Kiểm tra Tệp Mô tả Mẫu",
        f"Generated / Tạo lúc : {now}",
        f"Input / Đầu vào    : {args.input}",
        f"Status / Trạng thái: {status}",
        f"Samples / Số mẫu   : {len(validated_rows)}",
        "",
    ]
    if errors_total:
        report_lines.append("=== ERRORS / LỖI ===")
        report_lines.extend(errors_total)
        report_lines.append("")
    if warnings_total:
        report_lines.append("=== WARNINGS / CẢNH BÁO ===")
        report_lines.extend(warnings_total)
        report_lines.append("")
    if not errors_total and not warnings_total:
        report_lines.append("No issues found / Không có vấn đề nào.")

    with open(args.report, "w", encoding="utf-8") as fh:
        fh.write("\n".join(report_lines) + "\n")

    # Tóm tắt trạng thái / Summary status
    _progress(
        f"[validate_samplesheet] {status}: {len(validated_rows)} mẫu hợp lệ / valid sample(s); "
        f"{len(errors_total)} lỗi / error(s); {len(warnings_total)} cảnh báo / warning(s)"
    )
    for e in errors_total:
        print(e, file=sys.stderr)
    for w in warnings_total:
        print(w, file=sys.stderr)

    # ─── Bước 3: Ghi TSV đã kiểm tra / Step 3: Write validated TSV ─────────
    _progress(f"[validate_samplesheet] Ghi TSV đầu ra / Writing validated TSV: {args.output}")
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample_id", "ab1_forward", "ab1_reverse"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(validated_rows)

    _progress("[validate_samplesheet] Hoàn thành kiểm tra / Validation complete.")

    if errors_total:
        sys.exit(1)


if __name__ == "__main__":
    main()
