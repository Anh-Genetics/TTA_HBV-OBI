// ==============================================================================
// Module: ANNOTATE_VARIANTS
// [EN] Annotate variants in the per-sample consensus against the OBI-associated
//      mutation catalogue (data/references/obi_mutation_catalogue.tsv).
//
//      For each amino acid position in the S-ORF that differs from the
//      reference:
//        - Look up the position in the OBI catalogue.
//        - Report domain (a-determinant, MHR, S, pre-S2, pre-S1).
//        - Assign evidence level: A (strong), B (moderate), C (novel/weak).
//        - Flag premature stop codons and pre-S initiation codon changes.
//
//      NOTE: Domain coordinates are approximate (NC_003977.2, genotype A).
//            Alignment-based coordinate remapping is used for other genotypes.
//
// [VI] Chú thích biến thể trong trình tự đồng thuận mẫu dựa trên danh mục
//      đột biến liên quan OBI (data/references/obi_mutation_catalogue.tsv).
//
//      Với mỗi vị trí amino acid trong S-ORF khác với tham chiếu:
//        - Tra cứu vị trí trong danh mục OBI.
//        - Báo cáo vùng (a-determinant, MHR, S, pre-S2, pre-S1).
//        - Xếp loại bằng chứng: A (mạnh), B (vừa), C (mới/yếu).
//        - Gắn cờ codon dừng sớm và thay đổi codon khởi đầu pre-S.
//
//      LƯU Ý: Tọa độ vùng gần đúng (NC_003977.2, genotype A).
//              Ánh xạ tọa độ dựa trên alignment được dùng cho genotype khác.
// ==============================================================================

process ANNOTATE_VARIANTS {

    tag "${sample_id}"
    label 'process_medium'

    // [EN] Publish variant annotation TSV for report and manual review
    // [VI] Lưu TSV chú thích biến thể để báo cáo và xem xét thủ công
    publishDir "${params.outdir}/variants/${sample_id}", mode: 'copy'

    input:
    // [EN] Consensus FASTA from BUILD_CONSENSUS
    // [VI] FASTA đồng thuận từ BUILD_CONSENSUS
    tuple val(sample_id), path(consensus_fasta)

    // [EN] Multiple sequence alignment from ALIGN_CONSENSUS
    // [VI] Alignment đa trình tự từ ALIGN_CONSENSUS
    tuple val(sample_id), path(alignment_fasta)

    // [EN] OBI mutation catalogue TSV
    // [VI] Danh mục đột biến OBI dạng TSV
    path mutation_db

    output:
    // [EN] Per-variant annotation TSV (one row per variant)
    // [VI] TSV chú thích biến thể (một hàng mỗi biến thể)
    tuple val(sample_id), path("${sample_id}_variants.tsv"), emit: annotation_tsv

    script:
    """
    echo "[ANNOTATE_VARIANTS] ${sample_id}: Bắt đầu chú thích biến thể / Starting variant annotation..." >&2

    annotate_variants.py \\
        --sample-id   "${sample_id}" \\
        --consensus   "${consensus_fasta}" \\
        --alignment   "${alignment_fasta}" \\
        --mutation-db "${mutation_db}" \\
        --out-tsv     "${sample_id}_variants.tsv"

    # [EN] Count and report annotated variants
    # [VI] Đếm và báo cáo số biến thể được chú thích
    NVARS=\$(tail -n +2 "${sample_id}_variants.tsv" | wc -l)
    echo "[ANNOTATE_VARIANTS] ${sample_id}: \${NVARS} biến thể / variant(s) được chú thích / annotated." >&2
    """
}
