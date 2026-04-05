// ==============================================================================
// Module: ALIGN_CONSENSUS
// [EN] Align the per-sample consensus to the HBV reference panel using MAFFT.
//      The alignment is used downstream for:
//        - Domain-coordinate mapping in ANNOTATE_VARIANTS
//        - Manual phylogenetic analysis (run IQ-TREE separately on outputs)
//
// [VI] Căn chỉnh trình tự đồng thuận từng mẫu với bộ tham chiếu HBV dùng MAFFT.
//      Alignment được dùng tiếp theo cho:
//        - Ánh xạ tọa độ vùng trong ANNOTATE_VARIANTS
//        - Phân tích phân loài thủ công (chạy IQ-TREE riêng trên đầu ra)
//
// [EN] FALLBACK: If MAFFT is not installed, outputs a concatenated FASTA
//      (unaligned) with a warning – downstream annotation will use
//      approximate coordinates rather than alignment-derived positions.
// [VI] DỰ PHÒNG: Nếu MAFFT chưa cài đặt, xuất FASTA ghép (chưa căn chỉnh)
//      kèm cảnh báo – annotation sẽ dùng tọa độ gần đúng thay vì từ alignment.
// ==============================================================================

process ALIGN_CONSENSUS {

    tag "${sample_id}"
    label 'process_medium'

    // [EN] Publish alignment for manual inspection and phylogenetic analysis
    // [VI] Lưu alignment để kiểm tra thủ công và phân tích phân loài
    publishDir "${params.outdir}/alignments/${sample_id}", mode: 'copy'

    input:
    // [EN] Consensus FASTA from BUILD_CONSENSUS
    // [VI] FASTA đồng thuận từ BUILD_CONSENSUS
    tuple val(sample_id), path(consensus_fasta)

    // [EN] Shared HBV genotype reference FASTA
    // [VI] FASTA tham chiếu genotype HBV dùng chung
    path reference_fasta

    output:
    // [EN] Multiple sequence alignment in FASTA format
    // [VI] Alignment đa trình tự dạng FASTA
    tuple val(sample_id), path("${sample_id}_aligned.fasta"), emit: alignment_fasta

    script:
    """
    echo "[ALIGN_CONSENSUS] ${sample_id}: Bắt đầu căn chỉnh MAFFT / Starting MAFFT alignment..." >&2

    # [EN] Combine consensus + references for alignment input
    # [VI] Kết hợp đồng thuận + tham chiếu để tạo đầu vào alignment
    cat "${consensus_fasta}" "${reference_fasta}" > combined_for_alignment.fasta

    # [EN] Check if MAFFT is available; fall back to unaligned if not
    # [VI] Kiểm tra MAFFT có sẵn; dùng FASTA chưa căn chỉnh nếu không có
    if command -v mafft >/dev/null 2>&1; then
        echo "[ALIGN_CONSENSUS] ${sample_id}: Chạy MAFFT / Running MAFFT..." >&2
        mafft --auto --quiet combined_for_alignment.fasta > "${sample_id}_aligned.fasta"
        echo "[ALIGN_CONSENSUS] ${sample_id}: MAFFT hoàn thành / MAFFT complete." >&2
    else
        echo "[ALIGN_CONSENSUS] CẢNH BÁO / WARN: mafft không tìm thấy / not found." >&2
        echo "[ALIGN_CONSENSUS] Dùng FASTA chưa căn chỉnh / Using unaligned FASTA." >&2
        echo "[ALIGN_CONSENSUS] Cài đặt / Install: conda install -c bioconda mafft" >&2
        # [EN] Prepend a FASTA comment line so downstream processes can detect unaligned input
        # [VI] Thêm dòng ghi chú FASTA để bước sau nhận biết đầu vào chưa căn chỉnh
        echo "# WARNING: MAFFT not available – sequences below are UNALIGNED / CẢNH BÁO: chưa căn chỉnh" \\
            > "${sample_id}_aligned.fasta"
        cat combined_for_alignment.fasta >> "${sample_id}_aligned.fasta"
    fi
    """
}
