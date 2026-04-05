#!/usr/bin/env python3
"""
validate_samplesheet.py
=======================
Validate the HBV-OBI pipeline sample sheet CSV.

Expected CSV columns (case-insensitive, order-independent):
  sample_id     – unique sample identifier
  ab1_forward   – path to forward (.ab1) read
  ab1_reverse   – path to reverse (.ab1) read  [OPTIONAL; leave blank for single-read]

Outputs:
  --output  : validated TSV with canonical column names and resolved absolute paths
  --report  : human-readable validation report (text)

Exit codes:
  0  – validation passed (warnings may still be present)
  1  – one or more fatal errors found
"""

import argparse
import csv
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_COLUMNS = {"sample_id", "ab1_forward"}
OPTIONAL_COLUMNS = {"ab1_reverse"}


def parse_args():
    p = argparse.ArgumentParser(
        description="Validate the HBV-OBI pipeline sample sheet."
    )
    p.add_argument("--input",  required=True, help="Input CSV sample sheet")
    p.add_argument("--output", required=True, help="Output validated TSV")
    p.add_argument("--report", required=True, help="Output validation report text file")
    return p.parse_args()


def normalise_header(header: list[str]) -> dict[str, str]:
    """Return a mapping of lower-cased column names to original names."""
    return {col.strip().lower(): col.strip() for col in header}


def validate_file_path(path_str: str, field: str, sample_id: str) -> tuple[str, list[str]]:
    """
    Validate that a file path exists and has a recognised ABI extension.
    Returns (resolved_absolute_path, [error_messages]).
    """
    errors = []
    if not path_str or path_str.strip() == "":
        return "NONE", errors  # optional field is absent

    p = Path(path_str.strip())
    if not p.is_absolute():
        p = Path.cwd() / p  # resolve relative to CWD (Nextflow work dir)

    if not p.exists():
        errors.append(
            f"  [ERROR] {field} for sample '{sample_id}': file not found: {p}"
        )
    elif p.suffix.lower() not in {".ab1", ".abi"}:
        errors.append(
            f"  [WARN]  {field} for sample '{sample_id}': unexpected extension "
            f"'{p.suffix}' (expected .ab1 or .abi) – continuing anyway"
        )

    return str(p), errors


def main():
    args = parse_args()
    errors_total = []
    warnings_total = []
    validated_rows = []

    # -------------------------------------------------------------------------
    # Read input CSV
    # -------------------------------------------------------------------------
    try:
        with open(args.input, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None:
                sys.exit(f"[ERROR] Sample sheet '{args.input}' appears to be empty.")

            norm_map = normalise_header(list(reader.fieldnames))
            lower_cols = set(norm_map.keys())

            # Check required columns
            missing = REQUIRED_COLUMNS - lower_cols
            if missing:
                sys.exit(
                    f"[ERROR] Missing required column(s) in sample sheet: "
                    f"{', '.join(sorted(missing))}.\n"
                    f"Found columns: {', '.join(sorted(lower_cols))}"
                )

            has_reverse = "ab1_reverse" in lower_cols
            seen_ids: set[str] = set()

            for lineno, row in enumerate(reader, start=2):
                # Normalise keys
                norm_row = {k.strip().lower(): v for k, v in row.items()}

                sample_id = norm_row.get("sample_id", "").strip()
                if not sample_id:
                    errors_total.append(
                        f"  [ERROR] Line {lineno}: empty sample_id – skipping row"
                    )
                    continue

                if sample_id in seen_ids:
                    errors_total.append(
                        f"  [ERROR] Duplicate sample_id '{sample_id}' at line {lineno}"
                    )
                seen_ids.add(sample_id)

                # Forward read
                fwd_raw = norm_row.get("ab1_forward", "")
                fwd_path, fwd_errs = validate_file_path(fwd_raw, "ab1_forward", sample_id)
                if fwd_path == "NONE":
                    errors_total.append(
                        f"  [ERROR] Line {lineno}: sample '{sample_id}' has no ab1_forward value"
                    )
                else:
                    for e in fwd_errs:
                        if "[ERROR]" in e:
                            errors_total.append(e)
                        else:
                            warnings_total.append(e)

                # Reverse read (optional)
                rev_raw = norm_row.get("ab1_reverse", "") if has_reverse else ""
                rev_path, rev_errs = validate_file_path(rev_raw, "ab1_reverse", sample_id)
                for e in rev_errs:
                    if "[ERROR]" in e:
                        errors_total.append(e)
                    else:
                        warnings_total.append(e)

                validated_rows.append(
                    {
                        "sample_id": sample_id,
                        "ab1_forward": fwd_path,
                        "ab1_reverse": rev_path,  # "NONE" if absent
                    }
                )

    except FileNotFoundError:
        sys.exit(f"[ERROR] Sample sheet not found: {args.input}")
    except csv.Error as exc:
        sys.exit(f"[ERROR] Failed to parse CSV: {exc}")

    # -------------------------------------------------------------------------
    # Write validation report
    # -------------------------------------------------------------------------
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    status = "PASSED" if not errors_total else "FAILED"

    report_lines = [
        f"HBV-OBI Sample Sheet Validation Report",
        f"Generated: {now}",
        f"Input    : {args.input}",
        f"Status   : {status}",
        f"Samples  : {len(validated_rows)}",
        "",
    ]
    if errors_total:
        report_lines.append("=== ERRORS ===")
        report_lines.extend(errors_total)
        report_lines.append("")
    if warnings_total:
        report_lines.append("=== WARNINGS ===")
        report_lines.extend(warnings_total)
        report_lines.append("")
    if not errors_total and not warnings_total:
        report_lines.append("No issues found.")

    with open(args.report, "w", encoding="utf-8") as fh:
        fh.write("\n".join(report_lines) + "\n")

    # Print summary to stderr for Nextflow log visibility
    print(f"[validate_samplesheet] {status}: {len(validated_rows)} sample(s)", file=sys.stderr)
    for e in errors_total:
        print(e, file=sys.stderr)
    for w in warnings_total:
        print(w, file=sys.stderr)

    # -------------------------------------------------------------------------
    # Write validated TSV
    # -------------------------------------------------------------------------
    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample_id", "ab1_forward", "ab1_reverse"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(validated_rows)

    if errors_total:
        sys.exit(1)


if __name__ == "__main__":
    main()
