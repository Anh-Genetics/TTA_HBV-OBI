#!/usr/bin/env python3
"""
generate_report.py
==================
[EN] Generate a per-sample analysis report (HTML, TSV, or Markdown) summarising:
     - Read QC / trimming statistics
     - Consensus assembly quality
     - HBV genotype assignment
     - OBI-associated variant annotation with evidence classification
     - Clinical interpretation notes (for research use only)

[VI] Tạo báo cáo phân tích cho từng mẫu (HTML, TSV hoặc Markdown) tóm tắt:
     - Thống kê QC/cắt tỉa đoạn đọc
     - Chất lượng lắp ráp trình tự đồng thuận
     - Gán genotype HBV
     - Chú thích biến thể liên quan OBI kèm xếp hạng bằng chứng
     - Ghi chú diễn giải lâm sàng (chỉ dùng cho nghiên cứu)

[EN] Each sample runs in its own isolated process → no shared-file race conditions.
[VI] Mỗi mẫu chạy trong tiến trình riêng → không xảy ra xung đột ghi tệp.

[EN] Usage:
  generate_report.py --sample-id S001 --consensus S001_consensus.fasta \\
      --genotype S001_genotype.tsv --variants S001_variants.tsv \\
      --trim-stats S001_trim_stats.tsv --consensus-stats S001_consensus_stats.tsv \\
      --format html --out S001_report.html

[VI] Cách dùng:
  generate_report.py --sample-id M001 --consensus M001_consensus.fasta \\
      --genotype M001_genotype.tsv --variants M001_variants.tsv \\
      --trim-stats M001_trim_stats.tsv --consensus-stats M001_consensus_stats.tsv \\
      --format html --out M001_report.html
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
            "[EN] Generate per-sample HBV-OBI analysis report.\n"
            "[VI] Tạo báo cáo phân tích HBV-OBI cho từng mẫu."
        )
    )
    p.add_argument("--sample-id",       required=True, dest="sample_id",
                   help="[EN] Sample identifier / [VI] Mã mẫu")
    p.add_argument("--consensus",       required=True,
                   help="[EN] Consensus FASTA / [VI] FASTA đồng thuận")
    p.add_argument("--genotype",        required=True,
                   help="[EN] Genotype assignment TSV / [VI] TSV gán genotype")
    p.add_argument("--variants",        required=True,
                   help="[EN] Variant annotation TSV / [VI] TSV chú thích biến thể")
    p.add_argument("--trim-stats",      required=True, dest="trim_stats",
                   help="[EN] Trimming statistics TSV / [VI] TSV thống kê cắt tỉa")
    p.add_argument("--consensus-stats", required=True, dest="consensus_stats",
                   help="[EN] Consensus statistics TSV / [VI] TSV thống kê đồng thuận")
    p.add_argument("--format",          choices=["html", "tsv", "markdown"], default="html",
                   help="[EN] Output format / [VI] Định dạng đầu ra")
    p.add_argument("--out",             required=True,
                   help="[EN] Output report file / [VI] Tệp báo cáo đầu ra")
    return p.parse_args()


def read_tsv(path: str) -> list[dict]:
    """
    [EN] Read a TSV file into a list of dicts.
    [VI] Đọc tệp TSV thành danh sách dict.
    """
    rows = []
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows.append(row)
    except FileNotFoundError:
        pass
    return rows


def read_fasta_first(path: str) -> tuple[str, str]:
    """
    [EN] Return (header, sequence) of the first FASTA record.
    [VI] Trả về (tiêu_đề, trình_tự) của bản ghi FASTA đầu tiên.
    """
    header, parts = "", []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if not header:
                        header = line[1:]
                else:
                    parts.append(line)
    except FileNotFoundError:
        pass
    return header, "".join(parts)


def evidence_badge(level: str) -> str:
    """
    [EN] Return an HTML badge string for evidence level A/B/C.
         A=Strong/Mạnh, B=Moderate/Vừa, C=Weak-Novel/Yếu-Mới
    [VI] Trả về chuỗi HTML badge cho mức bằng chứng A/B/C.
    """
    colours = {"A": "#d32f2f", "B": "#f57c00", "C": "#388e3c"}
    labels  = {
        "A": "Strong / Mạnh",
        "B": "Moderate / Vừa",
        "C": "Weak/Novel / Yếu-Mới",
    }
    colour = colours.get(level.upper(), "#757575")
    label  = labels.get(level.upper(), level)
    return (
        f'<span style="background:{colour};color:#fff;padding:2px 6px;'
        f'border-radius:3px;font-size:0.85em;">{label}</span>'
    )


def generate_html(
    sample_id: str,
    cons_header: str,
    cons_seq: str,
    genotype_rows: list[dict],
    variant_rows: list[dict],
    trim_rows: list[dict],
    consensus_rows: list[dict],
    generated_at: str,
) -> str:
    """
    [EN] Generate a styled HTML report for a single sample.
    [VI] Tạo báo cáo HTML có định dạng cho một mẫu.
    """
    gt_row   = genotype_rows[0] if genotype_rows else {}
    gt       = gt_row.get("genotype", "N/A")
    subgt    = gt_row.get("subgenotype", "")
    gt_conf  = gt_row.get("confidence", "N/A")
    gt_pct   = gt_row.get("pct_identity", "N/A")

    cs_row       = consensus_rows[0] if consensus_rows else {}
    cons_len     = cs_row.get("consensus_length", len(cons_seq))
    cons_method  = cs_row.get("method", "N/A")
    n_conflicts  = cs_row.get("conflict_positions", "0")

    # ── Xây dựng bảng biến thể / Build variant table ─────────────────────────
    var_html = ""
    if variant_rows:
        for v in variant_rows:
            badge = evidence_badge(v.get("evidence_level", "C"))
            var_html += f"""
            <tr>
              <td>{v.get('notation','')}</td>
              <td>{v.get('domain','')}</td>
              <td>{v.get('mechanism','')}</td>
              <td>{badge}</td>
              <td>{v.get('catalogued','')}</td>
              <td>{v.get('notes','')}</td>
            </tr>"""
    else:
        var_html = "<tr><td colspan='6'><em>No variants annotated / Không có biến thể nào được chú thích</em></td></tr>"

    # ── Xây dựng bảng thống kê cắt tỉa / Build trim stats table ─────────────
    trim_html = ""
    for t in trim_rows:
        status_col = (
            '<td style="color:green">PASS ✅</td>'
            if t.get("qc_status") == "PASS"
            else '<td style="color:red">FAIL ❌</td>'
        )
        trim_html += f"""
        <tr>
          <td>{t.get('direction','')}</td>
          <td>{t.get('raw_length','')}</td>
          <td>{t.get('trimmed_length','')}</td>
          <td>{t.get('mean_quality','')}</td>
          <td>{t.get('ambiguous_bases','')}</td>
          {status_col}
        </tr>"""

    # ── Xác định cờ tổng quan / Determine overall flag ───────────────────────
    has_level_a = any(v.get("evidence_level") == "A" for v in variant_rows)
    overall_flag = (
        "⚠️ Cần xem xét / Review recommended"
        if has_level_a
        else "✅ Không phát hiện biến thể ưu tiên cao / No high-priority variants detected"
    )

    return f"""<!DOCTYPE html>
<html lang="vi">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>HBV-OBI Report – {sample_id}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2em; color: #212121; }}
    h1   {{ color: #1565C0; }}
    h2   {{ color: #1976D2; border-bottom: 1px solid #90CAF9; padding-bottom: 4px; }}
    table{{ border-collapse: collapse; width: 100%; margin-bottom: 1.5em; }}
    th   {{ background: #1565C0; color: #fff; padding: 8px 12px; text-align:left; }}
    td   {{ padding: 6px 12px; border-bottom: 1px solid #E0E0E0; }}
    tr:nth-child(even) {{ background: #F5F5F5; }}
    .flag-warn {{ background: #FFF3E0; border-left: 4px solid #F57C00;
                  padding: 12px; margin: 1em 0; }}
    .flag-ok   {{ background: #E8F5E9; border-left: 4px solid #388E3C;
                  padding: 12px; margin: 1em 0; }}
    pre  {{ background: #F5F5F5; padding: 12px; overflow-x: auto;
            font-size: 0.85em; white-space: pre-wrap; word-break: break-all; }}
    footer {{ color: #757575; font-size: 0.8em; margin-top: 3em; }}
    .disclaimer {{ background: #FFF8E1; border: 1px solid #FFD54F;
                   padding: 12px; font-size: 0.9em; margin-top: 2em; }}
  </style>
</head>
<body>

<h1>Báo cáo Phân tích HBV-OBI / HBV-OBI Sanger Analysis Report</h1>
<p><strong>Mã mẫu / Sample ID:</strong> {sample_id} &nbsp;|&nbsp;
   <strong>Tạo lúc / Generated:</strong> {generated_at} &nbsp;|&nbsp;
   <strong>Pipeline:</strong> TTA_HBV-OBI v0.1.0-MVP</p>

<div class="{'flag-warn' if has_level_a else 'flag-ok'}">
  {overall_flag}
</div>

<h2>1. Tóm tắt / Summary</h2>
<table>
  <tr><th>Thông số / Parameter</th><th>Giá trị / Value</th></tr>
  <tr><td>Genotype</td><td><strong>{gt}{(' / ' + subgt) if subgt else ''}</strong></td></tr>
  <tr><td>Độ tin cậy genotype / Genotype confidence</td><td>{gt_conf}</td></tr>
  <tr><td>% Đồng nhất BLAST / % Identity (BLAST)</td><td>{gt_pct}</td></tr>
  <tr><td>Độ dài đồng thuận / Consensus length (bp)</td><td>{cons_len}</td></tr>
  <tr><td>Phương pháp lắp ráp / Assembly method</td><td>{cons_method}</td></tr>
  <tr><td>Vị trí xung đột / Conflict positions</td><td>{n_conflicts}</td></tr>
  <tr><td>Biến thể được chú thích / Annotated variants</td><td>{len(variant_rows)}</td></tr>
</table>

<h2>2. QC / Cắt tỉa đoạn đọc / Read QC / Trimming</h2>
<table>
  <tr><th>Hướng / Direction</th><th>Độ dài thô / Raw (bp)</th>
      <th>Sau cắt / Trimmed (bp)</th><th>Q TB / Mean Q</th>
      <th>Base mơ hồ / Ambiguous</th><th>Trạng thái / Status</th></tr>
  {trim_html if trim_html else '<tr><td colspan="6"><em>Không có dữ liệu / No data</em></td></tr>'}
</table>

<h2>3. Gán Genotype / Genotype Assignment</h2>
<table>
  <tr><th>Genotype</th><th>Subgenotype</th><th>%ID</th>
      <th>BLAST hit tốt nhất / Best hit</th><th>Độ tin cậy / Confidence</th><th>Ghi chú / Note</th></tr>
  {''.join(
      f"<tr><td>{r.get('genotype','')}</td><td>{r.get('subgenotype','')}</td>"
      f"<td>{r.get('pct_identity','')}</td><td>{r.get('blast_hit','')}</td>"
      f"<td>{r.get('confidence','')}</td><td>{r.get('note','')}</td></tr>"
      for r in genotype_rows
  ) if genotype_rows else '<tr><td colspan="6"><em>Không có dữ liệu / No data</em></td></tr>'}
</table>

<h2>4. Chú thích biến thể OBI / OBI-Associated Variant Annotation</h2>
<table>
  <tr><th>Ký hiệu / Notation</th><th>Vùng / Domain</th>
      <th>Cơ chế / Mechanism</th><th>Bằng chứng / Evidence</th>
      <th>Trong danh mục / Catalogued</th><th>Ghi chú / Notes</th></tr>
  {var_html}
</table>

<h2>5. Trình tự đồng thuận / Consensus Sequence</h2>
<pre>{cons_seq if cons_seq else '(trống – QC đọc thất bại / empty – read QC failed)'}</pre>

<div class="disclaimer">
  <strong>⚠️ Tuyên bố miễn trách / Disclaimer &amp; Limitations (MVP)</strong><br>
  <strong>[VI]</strong> Báo cáo này được tạo bởi pipeline TTA_HBV-OBI v0.1.0-MVP
  chỉ dành cho mục đích nghiên cứu. Đây <strong>không phải</strong> công cụ chẩn đoán lâm sàng.<br>
  • Gán genotype dựa trên BLAST; không dựng cây phân loài tự động.<br>
  • Diễn giải biến thể OBI dựa trên tài liệu đã công bố; biến thể mới xếp mức C.<br>
  • Sanger sequencing không phát hiện được biến thể minor trong quasispecies.<br>
  • Tọa độ vùng HBV dựa trên NC_003977.2 (genotype A), gần đúng cho các genotype khác.<br><br>
  <strong>[EN]</strong> For research use only. Genotype is BLAST pre-screen only.
  OBI variants are literature-based; novel variants are Evidence Level C.
  Sanger cannot detect minor variants. Domain coordinates are approximate.
</div>

<footer>TTA_HBV-OBI pipeline &nbsp;|&nbsp; {generated_at}</footer>
</body>
</html>"""


def generate_tsv(
    sample_id: str,
    genotype_rows: list[dict],
    variant_rows: list[dict],
    trim_rows: list[dict],
    consensus_rows: list[dict],
    generated_at: str,
) -> str:
    """
    [EN] Generate a plain TSV report (sections separated by headers).
    [VI] Tạo báo cáo TSV đơn giản (các phần phân tách bằng tiêu đề).
    """
    lines = [
        f"# Báo cáo Phân tích HBV-OBI / TTA_HBV-OBI Analysis Report\t{sample_id}\t{generated_at}",
        "",
        "## GENOTYPE",
    ]
    if genotype_rows:
        lines.append("\t".join(genotype_rows[0].keys()))
        for r in genotype_rows:
            lines.append("\t".join(str(v) for v in r.values()))
    lines += ["", "## QC_TRIM_STATS / THỐNG KÊ CẮT TỈA"]
    if trim_rows:
        lines.append("\t".join(trim_rows[0].keys()))
        for r in trim_rows:
            lines.append("\t".join(str(v) for v in r.values()))
    lines += ["", "## CONSENSUS_STATS / THỐNG KÊ ĐỒNG THUẬN"]
    if consensus_rows:
        lines.append("\t".join(consensus_rows[0].keys()))
        for r in consensus_rows:
            lines.append("\t".join(str(v) for v in r.values()))
    lines += ["", "## VARIANTS / BIẾN THỂ"]
    if variant_rows:
        lines.append("\t".join(variant_rows[0].keys()))
        for r in variant_rows:
            lines.append("\t".join(str(v) for v in r.values()))
    else:
        lines.append("# Không có biến thể / No variants annotated")
    return "\n".join(lines) + "\n"


def generate_markdown(
    sample_id: str,
    cons_seq: str,
    genotype_rows: list[dict],
    variant_rows: list[dict],
    trim_rows: list[dict],
    consensus_rows: list[dict],
    generated_at: str,
) -> str:
    """
    [EN] Generate a Markdown report for a single sample.
    [VI] Tạo báo cáo Markdown cho một mẫu.
    """
    gt_row  = genotype_rows[0] if genotype_rows else {}
    gt      = gt_row.get("genotype", "N/A")
    subgt   = gt_row.get("subgenotype", "")
    gt_conf = gt_row.get("confidence", "N/A")

    def table_md(rows: list[dict]) -> str:
        """Tạo bảng Markdown / Build Markdown table."""
        if not rows:
            return "_Không có dữ liệu / No data_\n"
        keys      = list(rows[0].keys())
        header    = "| " + " | ".join(keys) + " |"
        separator = "| " + " | ".join("---" for _ in keys) + " |"
        body      = "\n".join(
            "| " + " | ".join(str(r.get(k, "")) for k in keys) + " |"
            for r in rows
        )
        return f"{header}\n{separator}\n{body}\n"

    return f"""# Báo cáo Phân tích HBV-OBI / HBV-OBI Sanger Analysis Report – {sample_id}

**Tạo lúc / Generated:** {generated_at}
**Pipeline:** TTA_HBV-OBI v0.1.0-MVP

---

## Tóm tắt / Summary

| Thông số / Parameter | Giá trị / Value |
|---|---|
| Genotype | **{gt}{' / ' + subgt if subgt else ''}** |
| Độ tin cậy / Confidence | {gt_conf} |
| Số biến thể / Variants | {len(variant_rows)} |
| Độ dài đồng thuận / Consensus length | {len(cons_seq)} bp |

---

## QC / Cắt tỉa đoạn đọc / Read Trimming

{table_md(trim_rows)}

## Gán genotype / Genotype Assignment

{table_md(genotype_rows)}

## Chú thích biến thể OBI / OBI Variant Annotation

{table_md(variant_rows) if variant_rows else '_Không có biến thể / No variants annotated._'}

## Trình tự đồng thuận / Consensus Sequence

```
{cons_seq if cons_seq else '(trống / empty)'}
```

---

> **Tuyên bố miễn trách / Disclaimer:** Chỉ dùng cho nghiên cứu / Research use only. See README for full limitations.
"""


def main():
    args         = parse_args()
    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    _progress(
        f"[generate_report] === Bắt đầu tạo báo cáo / Starting report generation: "
        f"{args.sample_id} (định dạng / format={args.format}) ==="
    )

    # ─── Bước 1: Đọc tất cả dữ liệu đầu vào / Step 1: Load all inputs ───────
    _progress(f"[generate_report] Đọc dữ liệu đầu vào / Loading inputs...")
    _, cons_seq      = read_fasta_first(args.consensus)
    genotype_rows    = read_tsv(args.genotype)
    variant_rows     = read_tsv(args.variants)
    trim_rows        = read_tsv(args.trim_stats)
    cons_rows        = read_tsv(args.consensus_stats)

    _progress(
        f"[generate_report] Dữ liệu: genotype={len(genotype_rows)} hàng, "
        f"biến thể / variants={len(variant_rows)}, "
        f"cắt tỉa / trim={len(trim_rows)} đoạn đọc / read(s)"
    )

    # Cảnh báo biến thể mức A / Warn about Level-A variants
    level_a = [v.get("notation","") for v in variant_rows if v.get("evidence_level") == "A"]
    if level_a:
        _progress(
            f"[generate_report] *** CẢNH BÁO / ALERT: {len(level_a)} biến thể mức A / "
            f"Level-A variant(s) detected: {', '.join(level_a)} ***"
        )

    # ─── Bước 2: Tạo nội dung báo cáo / Step 2: Generate report content ─────
    _progress(f"[generate_report] Tạo nội dung báo cáo / Generating {args.format} content...")
    fmt = args.format.lower()
    if fmt == "html":
        content = generate_html(
            args.sample_id, "", cons_seq,
            genotype_rows, variant_rows, trim_rows, cons_rows, generated_at,
        )
    elif fmt == "tsv":
        content = generate_tsv(
            args.sample_id, genotype_rows, variant_rows, trim_rows, cons_rows, generated_at,
        )
    else:  # markdown
        content = generate_markdown(
            args.sample_id, cons_seq,
            genotype_rows, variant_rows, trim_rows, cons_rows, generated_at,
        )

    # ─── Bước 3: Ghi tệp báo cáo / Step 3: Write report file ────────────────
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(content)

    _progress(
        f"[generate_report] === Báo cáo đã ghi / Report written: {args.out} "
        f"({fmt}, {len(content)} ký tự / chars) ==="
    )


if __name__ == "__main__":
    main()
