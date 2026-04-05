// ==============================================================================
// Module: BUILD_CONSENSUS
// [EN] Build a per-sample consensus sequence from all trimmed reads belonging
//      to the same sample_id.  Reads are grouped by sample_id (via groupTuple
//      in main.nf) before this process runs – preventing race conditions.
//
//      Single forward read → used directly as consensus.
//      Forward + Reverse   → RC the reverse, find overlap, merge with chosen
//                            conflict policy (iupac | majority | n).
//      No overlap found    → fall back to the longer PASS read.
//
// [VI] Xây dựng trình tự đồng thuận cho từng mẫu từ tất cả đoạn đọc đã cắt
//      thuộc cùng sample_id. Đoạn đọc được nhóm theo sample_id (qua groupTuple
//      trong main.nf) trước khi tiến trình này chạy – tránh xung đột song song.
//
//      Chỉ đọc xuôi     → dùng nguyên đoạn đọc làm đồng thuận.
//      Xuôi + Ngược     → RC đọc ngược, tìm vùng chồng, ghép với chính sách
//                         xử lý xung đột đã chọn (iupac | majority | n).
//      Không tìm thấy vùng chồng → dùng đoạn đọc dài hơn vượt QC.
// ==============================================================================

process BUILD_CONSENSUS {

    tag "${sample_id}"
    label 'process_medium'

    // [EN] Publish consensus + stats for downstream analysis and audit
    // [VI] Lưu đồng thuận + thống kê để phân tích tiếp theo và kiểm định
    publishDir "${params.outdir}/consensus/${sample_id}", mode: 'copy'

    input:
    // [EN] All trimmed FASTAs for this sample (grouped from TRIM_READS)
    // [VI] Tất cả FASTA đã cắt của mẫu này (được nhóm từ TRIM_READS)
    tuple val(sample_id), path(trimmed_fastas)

    output:
    // [EN] Consensus FASTA used by downstream alignment / annotation
    // [VI] FASTA đồng thuận dùng cho alignment / chú thích tiếp theo
    tuple val(sample_id), path("${sample_id}_consensus.fasta"),       emit: consensus_fasta

    // [EN] Consensus assembly statistics (method, overlap length, conflicts)
    // [VI] Thống kê lắp ráp đồng thuận (phương pháp, độ dài vùng chồng, xung đột)
    tuple val(sample_id), path("${sample_id}_consensus_stats.tsv"),   emit: consensus_stats

    script:
    """
    echo "[BUILD_CONSENSUS] ${sample_id}: Bắt đầu / Starting " \\
         "conflict=${params.conflict_policy} overlap=${params.min_overlap}" >&2

    build_consensus.py \\
        --sample-id   "${sample_id}" \\
        --reads       ${trimmed_fastas} \\
        --min-overlap "${params.min_overlap}" \\
        --conflict    "${params.conflict_policy}" \\
        --out-fasta   "${sample_id}_consensus.fasta" \\
        --out-stats   "${sample_id}_consensus_stats.tsv"

    echo "[BUILD_CONSENSUS] ${sample_id}: Hoàn thành / Done." >&2
    """
}
