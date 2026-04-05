#!/usr/bin/env python3
"""
annotate_variants.py
====================
[EN] Annotate variants in an HBV consensus sequence against:
     1. The OBI-associated mutation catalogue (evidence-graded A/B/C).
     2. HBV ORF/domain coordinate map (pre-S1, pre-S2, S, MHR, a-determinant).

[VI] Chú thích biến thể trong trình tự đồng thuận HBV dựa trên:
     1. Danh mục đột biến liên quan OBI (bằng chứng xếp hạng A/B/C).
     2. Bản đồ tọa độ ORF/vùng HBV (pre-S1, pre-S2, S, MHR, vùng a-determinant).

[EN] Steps:
     1. Read MAFFT alignment (consensus + references).
     2. Identify reference row (by accession patterns).
     3. Translate S, pre-S2, pre-S1 ORFs in the consensus.
     4. Identify AA changes vs reference.
     5. Look up in OBI catalogue; output one TSV row per variant.

[VI] Các bước:
     1. Đọc alignment MAFFT (đồng thuận + tham chiếu).
     2. Xác định hàng tham chiếu (theo mẫu số hiệu).
     3. Dịch mã ORF S, pre-S2, pre-S1 trong đồng thuận.
     4. Xác định thay đổi amino acid so với tham chiếu.
     5. Tra cứu trong danh mục OBI; xuất một hàng TSV mỗi biến thể.

[EN] NOTE (MVP): Domain coordinates are approximate (genotype A, NC_003977.2).
[VI] LƯU Ý (MVP): Tọa độ vùng là gần đúng (genotype A, NC_003977.2).

[EN] Usage:
  annotate_variants.py --sample-id S001 --consensus S001_consensus.fasta \\
    --alignment S001_aligned.fasta --mutation-db catalogue.tsv \\
    --out-tsv S001_variants.tsv

[VI] Cách dùng:
  annotate_variants.py --sample-id M001 --consensus M001_consensus.fasta \\
    --alignment M001_aligned.fasta --mutation-db catalogue.tsv \\
    --out-tsv M001_variants.tsv
"""

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

# Kiểm tra Biopython / Check Biopython
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "[ERROR] Biopython chưa được cài đặt / Biopython not found.\n"
        "Cài đặt / Install: conda install -c conda-forge biopython"
    )

# ─── Bản đồ tọa độ vùng HBV (tính gần đúng, dựa trên amplicon ~1245 bp)
# [EN] HBV domain coordinate map (approximate, based on ~1245 bp preS1-S amplicon)
# [VI] Bản đồ tọa độ vùng HBV (gần đúng, dựa trên amplicon preS1-S ~1245 bp)
# Định dạng: (vị_trí_bắt_đầu_nt_0_based, vị_trí_kết_thúc_exclusive)
DOMAIN_MAP = {
    "preS1":         (0,    400),
    "preS2":         (400,  550),
    "S":             (550, 1245),
    "MHR":           (850, 1057),   # S codon ~99-169
    "a_determinant": (922,  991),   # S codon ~124-147
}


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
            "[EN] Annotate HBV consensus variants against OBI catalogue.\n"
            "[VI] Chú thích biến thể đồng thuận HBV theo danh mục OBI."
        )
    )
    p.add_argument("--sample-id",   required=True, dest="sample_id")
    p.add_argument("--consensus",   required=True,
                   help="[EN] Per-sample consensus FASTA / [VI] FASTA đồng thuận mẫu")
    p.add_argument("--alignment",   required=True,
                   help="[EN] MAFFT alignment FASTA (consensus+refs) / [VI] Alignment MAFFT")
    p.add_argument("--mutation-db", required=True, dest="mutation_db",
                   help="[EN] OBI mutation catalogue TSV / [VI] Danh mục đột biến OBI TSV")
    p.add_argument("--out-tsv",     required=True, dest="out_tsv",
                   help="[EN] Output variant annotation TSV / [VI] TSV chú thích biến thể đầu ra")
    return p.parse_args()


def read_fasta_dict(path: str) -> dict[str, str]:
    """
    [EN] Return {seq_id: sequence_no_gaps} from a FASTA file.
    [VI] Trả về {seq_id: trình_tự_không_gap} từ tệp FASTA.
    """
    records = {}
    try:
        for rec in SeqIO.parse(path, "fasta"):
            records[rec.id] = str(rec.seq).upper()
    except Exception as exc:  # noqa: BLE001
        _progress(f"[annotate_variants] CẢNH BÁO / WARN: Không đọc được FASTA '{path}': {exc}")
    return records


def read_alignment(path: str) -> dict[str, str]:
    """
    [EN] Return {seq_id: gapped_sequence} from a multiple alignment FASTA.
    [VI] Trả về {seq_id: trình_tự_có_gap} từ tệp alignment FASTA.
    """
    records = {}
    try:
        for rec in SeqIO.parse(path, "fasta"):
            records[rec.id] = str(rec.seq).upper()
    except Exception as exc:  # noqa: BLE001
        _progress(f"[annotate_variants] CẢNH BÁO / WARN: Không đọc được alignment '{path}': {exc}")
    return records


def load_mutation_catalogue(path: str) -> list[dict]:
    """
    [EN] Load the OBI mutation catalogue TSV.
         Expected columns: aa_position, ref_aa, alt_aa, region, mechanism,
                           evidence_level, notes
    [VI] Tải danh mục đột biến OBI từ tệp TSV.
         Cột mong đợi: aa_position, ref_aa, alt_aa, region, mechanism,
                        evidence_level, notes
    """
    catalogue = []
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    row["aa_position"] = int(row["aa_position"])
                except (KeyError, ValueError):
                    pass
                catalogue.append(row)
    except FileNotFoundError:
        _progress(f"[annotate_variants] CẢNH BÁO / WARN: Không tìm thấy danh mục / Catalogue not found: {path}")
    return catalogue


def find_reference_id(alignment: dict[str, str], sample_id: str) -> str | None:
    """
    [EN] Pick the reference sequence from the alignment (not the sample).
         Prefers sequences with accession-like IDs (NC_, AB_, AY_, JN_...).
    [VI] Chọn trình tự tham chiếu từ alignment (không phải mẫu).
         Ưu tiên ID dạng số hiệu (NC_, AB_, AY_, JN_...).
    """
    for seq_id in alignment:
        if seq_id == sample_id:
            continue
        if any(kw in seq_id.upper() for kw in ("REF", "REFERENCE", "NC_", "AB_", "AY_", "JN_")):
            return seq_id
    # Dự phòng: chuỗi không phải mẫu đầu tiên / Fallback: first non-sample sequence
    for seq_id in alignment:
        if seq_id != sample_id:
            return seq_id
    return None


def translate_frame(nt_seq: str, frame: int = 0) -> str:
    """
    [EN] Translate a nucleotide string starting at *frame* offset.
    [VI] Dịch mã chuỗi nucleotide bắt đầu từ offset *frame*.
    """
    sub = nt_seq[frame:]
    remainder = len(sub) % 3
    if remainder:
        sub += "N" * (3 - remainder)  # Đệm để chia hết 3 / Pad to multiple of 3
    try:
        return str(Seq(sub).translate(to_stop=False))
    except Exception:  # noqa: BLE001
        return ""


def _domain_label(nt_pos_in_amplicon: int) -> str:
    """
    [EN] Return the most specific domain label for a nucleotide position.
         Priority order: a_determinant > MHR > S > preS2 > preS1
    [VI] Trả về nhãn vùng cụ thể nhất cho một vị trí nucleotide.
         Thứ tự ưu tiên: a_determinant > MHR > S > preS2 > preS1
    """
    priority_order = ["a_determinant", "MHR", "S", "preS2", "preS1"]
    for name in priority_order:
        start, end = DOMAIN_MAP[name]
        if start <= nt_pos_in_amplicon < end:
            return name
    return "outside_amplicon"


def annotate(
    sample_id: str,
    consensus_seq: str,
    ref_seq: str,
    catalogue: list[dict],
) -> list[dict]:
    """
    [EN] Compare consensus_seq vs ref_seq (both ungapped) and look up variants
         in the OBI catalogue. Returns list of variant annotation dicts.
    [VI] So sánh consensus_seq với ref_seq (cả hai không có gap) và tra cứu
         biến thể trong danh mục OBI. Trả về danh sách dict chú thích biến thể.
    """
    variants = []

    # ── Dịch mã vùng S / Translate S region ──────────────────────────────────
    s_start  = DOMAIN_MAP["S"][0]
    cons_s   = consensus_seq[s_start:] if len(consensus_seq) > s_start else ""
    ref_s    = ref_seq[s_start:]       if len(ref_seq) > s_start else ""
    cons_aa  = translate_frame(cons_s, 0) if cons_s else ""
    ref_aa   = translate_frame(ref_s,  0) if ref_s  else ""

    # ── Xây dựng tra cứu theo vị trí / Build catalogue lookup by position ────
    cat_by_pos: dict[int, list[dict]] = {}
    for entry in catalogue:
        pos = entry.get("aa_position")
        if isinstance(pos, int):
            cat_by_pos.setdefault(pos, []).append(entry)

    # ── Quét thay đổi amino acid / Scan amino acid changes ───────────────────
    max_len = min(len(cons_aa), len(ref_aa))
    for aa_pos in range(max_len):
        c_aa = cons_aa[aa_pos]
        r_aa = ref_aa[aa_pos]
        if c_aa == r_aa:
            continue

        pos_1based = aa_pos + 1
        domain     = _domain_label(s_start + aa_pos * 3)

        # Tra cứu trong danh mục / Lookup in catalogue
        cat_entries = cat_by_pos.get(pos_1based, [])
        matched     = [e for e in cat_entries if e.get("alt_aa", "").upper() == c_aa]

        if matched:
            for e in matched:
                variants.append(
                    {
                        "sample_id":      sample_id,
                        "aa_position":    pos_1based,
                        "ref_aa":         r_aa,
                        "alt_aa":         c_aa,
                        "notation":       f"s{r_aa}{pos_1based}{c_aa}",
                        "domain":         domain,
                        "mechanism":      e.get("mechanism", ""),
                        "evidence_level": e.get("evidence_level", ""),
                        "notes":          e.get("notes", ""),
                        "catalogued":     "yes",
                    }
                )
        else:
            # Biến thể mới / không có trong danh mục / Novel / uncatalogued variant
            variants.append(
                {
                    "sample_id":      sample_id,
                    "aa_position":    pos_1based,
                    "ref_aa":         r_aa,
                    "alt_aa":         c_aa,
                    "notation":       f"s{r_aa}{pos_1based}{c_aa}",
                    "domain":         domain,
                    "mechanism":      "unknown",
                    "evidence_level": "C",   # Mới / Novel
                    "notes":          "not_in_obi_catalogue",
                    "catalogued":     "no",
                }
            )

    # ── Kiểm tra codon dừng sớm / Check for premature stop codons ─────────────
    stop_positions    = [i + 1 for i, aa in enumerate(cons_aa) if aa == "*"]
    expected_s_len    = len(ref_aa)
    premature_stops   = [p for p in stop_positions if p < expected_s_len]

    for p in premature_stops:
        variants.append(
            {
                "sample_id":      sample_id,
                "aa_position":    p,
                "ref_aa":         ref_aa[p - 1] if p <= len(ref_aa) else "?",
                "alt_aa":         "*",
                "notation":       f"s{ref_aa[p-1] if p <= len(ref_aa) else '?'}{p}*",
                "domain":         _domain_label(s_start + (p - 1) * 3),
                "mechanism":      "truncation_of_HBsAg",
                "evidence_level": "A",
                "notes":          "premature_stop_codon_in_S_ORF",
                "catalogued":     "yes",
            }
        )

    return variants


def main():
    args = parse_args()

    _progress(
        f"[annotate_variants] === Bắt đầu chú thích biến thể / Starting variant annotation: "
        f"{args.sample_id} ==="
    )

    # ─── Bước 1: Đọc dữ liệu vào / Step 1: Load inputs ──────────────────────
    _progress(f"[annotate_variants] Đọc trình tự đồng thuận / Reading consensus: {args.consensus}")
    consensus_seqs = read_fasta_dict(args.consensus)

    _progress(f"[annotate_variants] Đọc alignment / Reading alignment: {args.alignment}")
    alignment = read_alignment(args.alignment)

    _progress(f"[annotate_variants] Đọc danh mục đột biến / Loading catalogue: {args.mutation_db}")
    catalogue = load_mutation_catalogue(args.mutation_db)
    _progress(f"[annotate_variants] Danh mục / Catalogue: {len(catalogue)} mục nhập / entries")

    # ─── Bước 2: Xác định trình tự đồng thuận / Step 2: Identify consensus ───
    cons_seq = ""
    for sid, seq in consensus_seqs.items():
        if args.sample_id in sid or sid == args.sample_id:
            cons_seq = seq.replace("-", "")
            break
    if not cons_seq and consensus_seqs:
        cons_seq = list(consensus_seqs.values())[0].replace("-", "")
    _progress(f"[annotate_variants] Đồng thuận / Consensus: {len(cons_seq)} bp")

    # ─── Bước 3: Xác định tham chiếu / Step 3: Identify reference ────────────
    ref_id  = find_reference_id(alignment, args.sample_id)
    ref_seq = alignment.get(ref_id, "").replace("-", "") if ref_id else ""
    _progress(
        f"[annotate_variants] Tham chiếu / Reference: "
        f"'{ref_id}' ({len(ref_seq)} bp)"
    )

    # ─── Bước 4: Chú thích biến thể / Step 4: Annotate variants ─────────────
    if not cons_seq:
        _progress(
            f"[annotate_variants] CẢNH BÁO / WARN: {args.sample_id}: "
            "trình tự đồng thuận rỗng – bỏ qua chú thích / "
            "empty consensus – skipping annotation"
        )
        variants = []
    elif not ref_seq:
        _progress(
            f"[annotate_variants] CẢNH BÁO / WARN: {args.sample_id}: "
            "không tìm thấy tham chiếu trong alignment – bỏ qua / "
            "no reference in alignment – skipping"
        )
        variants = []
    else:
        _progress("[annotate_variants] Dịch mã và so sánh amino acid / Translating and comparing amino acids...")
        variants = annotate(args.sample_id, cons_seq, ref_seq, catalogue)

    # ─── Bước 5: Ghi TSV / Step 5: Write TSV ────────────────────────────────
    fieldnames = [
        "sample_id", "aa_position", "ref_aa", "alt_aa", "notation",
        "domain", "mechanism", "evidence_level", "notes", "catalogued",
    ]
    with open(args.out_tsv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(variants)

    # Tóm tắt kết quả / Result summary
    n_a = sum(1 for v in variants if v.get("evidence_level") == "A")
    n_b = sum(1 for v in variants if v.get("evidence_level") == "B")
    n_c = sum(1 for v in variants if v.get("evidence_level") == "C")

    _progress(
        f"[annotate_variants] Kết quả / Results: {len(variants)} biến thể / variant(s) | "
        f"Mức A={n_a}, B={n_b}, C={n_c}"
    )
    if n_a > 0:
        a_variants = [v["notation"] for v in variants if v.get("evidence_level") == "A"]
        _progress(
            f"[annotate_variants] *** CẢNH BÁO / REVIEW REQUIRED: "
            f"Biến thể mức A / Level-A variants: {', '.join(a_variants)} ***"
        )

    _progress(
        f"[annotate_variants] === Hoàn thành chú thích / Annotation complete: "
        f"{args.sample_id} ==="
    )


if __name__ == "__main__":
    main()
