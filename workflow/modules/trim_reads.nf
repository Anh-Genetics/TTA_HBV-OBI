// ==============================================================================
// Module: TRIM_READS
// [EN] Quality-trim a parsed Sanger read using its companion quality TSV.
//      Strategy: fixed-end trim + sliding-window quality scan from both ends.
//      Reads shorter than min_length after trimming are flagged FAIL and
//      produce an empty FASTA (downstream processes handle gracefully).
//
// [VI] Cắt tỉa chất lượng đoạn đọc Sanger đã phân tích dùng TSV chất lượng.
//      Chiến lược: cắt đầu cố định + quét chất lượng từ hai đầu.
//      Đoạn đọc quá ngắn sau cắt bị đánh dấu FAIL và tạo FASTA rỗng
//      (các bước sau xử lý gracefully).
// ==============================================================================

process TRIM_READS {

    tag "${sample_id}:${direction}"
    label 'process_low'

    // [EN] Save trimmed FASTA and stats for QC review
    // [VI] Lưu FASTA đã cắt và thống kê để kiểm tra QC
    publishDir "${params.outdir}/reads/trimmed/${sample_id}", mode: 'copy'

    input:
    // [EN] FASTA from PARSE_ABI
    // [VI] FASTA từ PARSE_ABI
    tuple val(sample_id), val(direction), path(raw_fasta)
    // [EN] Quality TSV from PARSE_ABI (joined by sample_id + direction)
    // [VI] TSV chất lượng từ PARSE_ABI (ghép theo mã mẫu và hướng)
    tuple val(sample_id), val(direction), path(quality_tsv)

    output:
    // [EN] Trimmed FASTA + trim statistics, all keyed by sample_id
    // [VI] FASTA đã cắt + thống kê cắt tỉa, theo mã mẫu
    tuple val(sample_id), val(direction), path("${sample_id}_${direction}_trimmed.fasta"),
          path("${sample_id}_${direction}_trim_stats.tsv"),  emit: trimmed_fasta
    // [EN] Collected trim stats for per-sample summary
    // [VI] Thống kê cắt tỉa tổng hợp cho tóm tắt mẫu
    tuple val(sample_id), path("${sample_id}_${direction}_trim_stats.tsv"),  emit: trim_summary

    script:
    """
    echo "[TRIM_READS] ${sample_id} (${direction}): Cắt tỉa / Trimming " \\
         "min_quality=${params.min_quality} min_length=${params.min_length}" >&2

    trim_reads.py \\
        --fasta       "${raw_fasta}" \\
        --quality     "${quality_tsv}" \\
        --sample-id   "${sample_id}" \\
        --direction   "${direction}" \\
        --min-quality "${params.min_quality}" \\
        --min-length  "${params.min_length}" \\
        --trim-ends   "${params.trim_ends_bases}" \\
        --out-fasta   "${sample_id}_${direction}_trimmed.fasta" \\
        --out-stats   "${sample_id}_${direction}_trim_stats.tsv"

    echo "[TRIM_READS] ${sample_id} (${direction}): Hoàn thành cắt tỉa / Trimming complete." >&2
    """
}
