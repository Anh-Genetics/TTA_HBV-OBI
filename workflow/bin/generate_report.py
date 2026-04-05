#!/usr/bin/env python3
"""
generate_report.py
==================
Generate a per-sample analysis report in HTML, TSV, or Markdown format.

Inputs
------
  --consensus      : per-sample consensus FASTA
  --genotype       : genotype assignment TSV
  --variants       : variant annotation TSV
  --trim-stats     : trimming statistics TSV
  --consensus-stats: consensus building statistics TSV
  --format         : output format (html | tsv | markdown)

Output
------
  A single file per sample (HTML, TSV, or Markdown) summarising all analysis
  results in a structured, reviewable format.

Usage:
  generate_report.py \\
      --sample-id S001 \\
      --consensus S001_consensus.fasta \\
      --genotype S001_genotype.tsv \\
      --variants S001_variants.tsv \\
      --trim-stats S001_trim_stats.tsv \\
      --consensus-stats S001_consensus_stats.tsv \\
      --format html \\
      --out S001_report.html
"""

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Generate per-sample HBV-OBI analysis report.")
    p.add_argument("--sample-id",       required=True, dest="sample_id")
    p.add_argument("--consensus",       required=True)
    p.add_argument("--genotype",        required=True)
    p.add_argument("--variants",        required=True)
    p.add_argument("--trim-stats",      required=True, dest="trim_stats")
    p.add_argument("--consensus-stats", required=True, dest="consensus_stats")
    p.add_argument("--format",          choices=["html", "tsv", "markdown"], default="html")
    p.add_argument("--out",             required=True)
    return p.parse_args()


def read_tsv(path: str) -> list[dict]:
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
    """Return (header, sequence) of the first record."""
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
    """Return an HTML badge string for evidence level A/B/C."""
    colours = {"A": "#d32f2f", "B": "#f57c00", "C": "#388e3c"}
    labels = {"A": "Strong", "B": "Moderate", "C": "Weak/Novel"}
    colour = colours.get(level.upper(), "#757575")
    label = labels.get(level.upper(), level)
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
    gt_row = genotype_rows[0] if genotype_rows else {}
    gt = gt_row.get("genotype", "N/A")
    subgt = gt_row.get("subgenotype", "")
    gt_conf = gt_row.get("confidence", "N/A")
    gt_pct = gt_row.get("pct_identity", "N/A")

    cs_row = consensus_rows[0] if consensus_rows else {}
    cons_len = cs_row.get("consensus_length", len(cons_seq))
    cons_method = cs_row.get("method", "N/A")
    n_conflicts = cs_row.get("conflict_positions", "0")

    # Build variant table rows
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
        var_html = "<tr><td colspan='6'><em>No variants annotated</em></td></tr>"

    # Build trim stats
    trim_html = ""
    for t in trim_rows:
        status_col = (
            '<td style="color:green">PASS</td>'
            if t.get("qc_status") == "PASS"
            else '<td style="color:red">FAIL</td>'
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

    overall_flag = (
        "⚠️ Review recommended"
        if any(v.get("evidence_level") == "A" for v in variant_rows)
        else "✅ No high-priority variants detected"
    )

    return f"""<!DOCTYPE html>
<html lang="en">
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

<h1>HBV-OBI Sanger Analysis Report</h1>
<p><strong>Sample ID:</strong> {sample_id} &nbsp;|&nbsp;
   <strong>Generated:</strong> {generated_at} &nbsp;|&nbsp;
   <strong>Pipeline:</strong> TTA_HBV-OBI v0.1.0-MVP</p>

<div class="{'flag-warn' if '⚠️' in overall_flag else 'flag-ok'}">
  {overall_flag}
</div>

<h2>1. Summary</h2>
<table>
  <tr><th>Parameter</th><th>Value</th></tr>
  <tr><td>Genotype</td><td><strong>{gt}{(' / ' + subgt) if subgt else ''}</strong></td></tr>
  <tr><td>Genotype confidence</td><td>{gt_conf}</td></tr>
  <tr><td>% Identity (BLAST)</td><td>{gt_pct}</td></tr>
  <tr><td>Consensus length (bp)</td><td>{cons_len}</td></tr>
  <tr><td>Assembly method</td><td>{cons_method}</td></tr>
  <tr><td>Conflict positions</td><td>{n_conflicts}</td></tr>
  <tr><td>Annotated variants</td><td>{len(variant_rows)}</td></tr>
</table>

<h2>2. Read QC / Trimming</h2>
<table>
  <tr><th>Direction</th><th>Raw length</th><th>Trimmed length</th>
      <th>Mean Q</th><th>Ambiguous bases</th><th>Status</th></tr>
  {trim_html if trim_html else '<tr><td colspan="6"><em>No trim stats available</em></td></tr>'}
</table>

<h2>3. Genotype Assignment</h2>
<table>
  <tr><th>Genotype</th><th>Subgenotype</th><th>% Identity</th>
      <th>Best BLAST hit</th><th>Confidence</th><th>Note</th></tr>
  {''.join(
      f"<tr><td>{r.get('genotype','')}</td><td>{r.get('subgenotype','')}</td>"
      f"<td>{r.get('pct_identity','')}</td><td>{r.get('blast_hit','')}</td>"
      f"<td>{r.get('confidence','')}</td><td>{r.get('note','')}</td></tr>"
      for r in genotype_rows
  ) if genotype_rows else '<tr><td colspan="6"><em>No genotype data</em></td></tr>'}
</table>

<h2>4. OBI-Associated Variant Annotation</h2>
<table>
  <tr><th>Notation</th><th>Domain</th><th>Mechanism</th>
      <th>Evidence</th><th>Catalogued</th><th>Notes</th></tr>
  {var_html}
</table>

<h2>5. Consensus Sequence</h2>
<pre>{cons_seq if cons_seq else '(empty – read QC failed)'}</pre>

<div class="disclaimer">
  <strong>⚠️ Disclaimer / Limitations (MVP)</strong><br>
  This report is generated by the TTA_HBV-OBI pipeline v0.1.0-MVP and is intended
  for research purposes only. It is <strong>not a clinical diagnostic tool</strong>.<br>
  • Genotype assignment is based on BLAST alignment against a curated reference panel;
    phylogenetic tree confirmation is not performed automatically.<br>
  • OBI variant interpretation is based on published literature; novel variants are
    flagged as Evidence Level C (weak/novel) and require independent validation.<br>
  • Minor variants and quasispecies are not detectable from Sanger data alone.<br>
  • HBV domain coordinates are approximate (based on genotype A reference NC_003977.2).
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
    lines = [
        f"# TTA_HBV-OBI Analysis Report\t{sample_id}\t{generated_at}",
        "",
        "## GENOTYPE",
    ]
    if genotype_rows:
        hdr = "\t".join(genotype_rows[0].keys())
        lines.append(hdr)
        for r in genotype_rows:
            lines.append("\t".join(str(v) for v in r.values()))
    lines += ["", "## QC_TRIM_STATS"]
    if trim_rows:
        lines.append("\t".join(trim_rows[0].keys()))
        for r in trim_rows:
            lines.append("\t".join(str(v) for v in r.values()))
    lines += ["", "## CONSENSUS_STATS"]
    if consensus_rows:
        lines.append("\t".join(consensus_rows[0].keys()))
        for r in consensus_rows:
            lines.append("\t".join(str(v) for v in r.values()))
    lines += ["", "## VARIANTS"]
    if variant_rows:
        lines.append("\t".join(variant_rows[0].keys()))
        for r in variant_rows:
            lines.append("\t".join(str(v) for v in r.values()))
    else:
        lines.append("# No variants annotated")
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
    gt_row = genotype_rows[0] if genotype_rows else {}
    gt = gt_row.get("genotype", "N/A")
    subgt = gt_row.get("subgenotype", "")
    gt_conf = gt_row.get("confidence", "N/A")

    def table_md(rows: list[dict]) -> str:
        if not rows:
            return "_No data_\n"
        keys = list(rows[0].keys())
        header = "| " + " | ".join(keys) + " |"
        separator = "| " + " | ".join("---" for _ in keys) + " |"
        body = "\n".join(
            "| " + " | ".join(str(r.get(k, "")) for k in keys) + " |"
            for r in rows
        )
        return f"{header}\n{separator}\n{body}\n"

    return f"""# HBV-OBI Sanger Analysis Report – {sample_id}

**Generated:** {generated_at}  
**Pipeline:** TTA_HBV-OBI v0.1.0-MVP

---

## Summary

| Parameter | Value |
|---|---|
| Genotype | **{gt}{' / ' + subgt if subgt else ''}** |
| Confidence | {gt_conf} |
| Variants | {len(variant_rows)} |
| Consensus length | {len(cons_seq)} bp |

---

## Read QC / Trimming

{table_md(trim_rows)}

## Genotype Assignment

{table_md(genotype_rows)}

## OBI Variant Annotation

{table_md(variant_rows) if variant_rows else '_No variants annotated._'}

## Consensus Sequence

```
{cons_seq if cons_seq else '(empty)'}
```

---

> **Disclaimer:** Research use only. See README for full limitations.
"""


def main():
    args = parse_args()
    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    _, cons_seq = read_fasta_first(args.consensus)
    genotype_rows = read_tsv(args.genotype)
    variant_rows  = read_tsv(args.variants)
    trim_rows     = read_tsv(args.trim_stats)
    cons_rows     = read_tsv(args.consensus_stats)

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

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(content)

    print(
        f"[generate_report] {args.sample_id}: report written to {args.out} ({fmt})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
