// ==============================================================================
// Module: GENERATE_REPORT
// [EN] Generate a per-sample analysis report summarising all pipeline outputs.
//      Output format is controlled by params.report_format (html|tsv|markdown).
//
//      IMPORTANT: Each sample runs in its own isolated Nextflow work directory.
//      There is NO shared file being written concurrently – no race conditions.
//
// [VI] Tạo báo cáo phân tích cho từng mẫu, tóm tắt tất cả đầu ra pipeline.
//      Định dạng đầu ra được kiểm soát bởi params.report_format (html|tsv|markdown).
//
//      QUAN TRỌNG: Mỗi mẫu chạy trong thư mục làm việc Nextflow riêng biệt.
//      KHÔNG có tệp dùng chung nào được ghi đồng thời – không xung đột.
//
// [EN] Inputs joined by sample_id in main.nf before reaching this process.
// [VI] Đầu vào được join theo sample_id trong main.nf trước khi đến tiến trình này.
// ==============================================================================

process GENERATE_REPORT {

    tag "${sample_id}"
    label 'process_low'

    // [EN] Save reports to a dedicated results directory by format
    // [VI] Lưu báo cáo vào thư mục kết quả riêng theo định dạng
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    // [EN] All per-sample outputs joined by sample_id:
    //      (sample_id, consensus_fasta, genotype_tsv, variant_tsv,
    //       trim_stats_tsv, consensus_stats_tsv)
    // [VI] Tất cả đầu ra từng mẫu join theo sample_id:
    //      (mã_mẫu, fasta_đồng_thuận, tsv_genotype, tsv_biến_thể,
    //       tsv_thống_kê_cắt, tsv_thống_kê_đồng_thuận)
    tuple val(sample_id),
          path(consensus_fasta),
          path(genotype_tsv),
          path(variant_tsv),
          path(trim_stats_tsv),
          path(consensus_stats_tsv)

    output:
    // [EN] Individual sample report file (HTML / TSV / Markdown)
    // [VI] Tệp báo cáo mẫu đơn lẻ (HTML / TSV / Markdown)
    path "${sample_id}_report.${params.report_format}", emit: sample_report

    script:
    """
    echo "[GENERATE_REPORT] ${sample_id}: Tạo báo cáo / Generating report (${params.report_format})..." >&2

    generate_report.py \\
        --sample-id       "${sample_id}" \\
        --consensus       "${consensus_fasta}" \\
        --genotype        "${genotype_tsv}" \\
        --variants        "${variant_tsv}" \\
        --trim-stats      "${trim_stats_tsv}" \\
        --consensus-stats "${consensus_stats_tsv}" \\
        --format          "${params.report_format}" \\
        --out             "${sample_id}_report.${params.report_format}"

    echo "[GENERATE_REPORT] ${sample_id}: Báo cáo đã ghi / Report written: " \\
         "${sample_id}_report.${params.report_format}" >&2
    """
}
