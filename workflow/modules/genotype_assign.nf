// ==============================================================================
// Module: ASSIGN_GENOTYPE
// [EN] Assign HBV genotype (A–I) to the per-sample consensus using BLAST.
//      The consensus is BLASTed against the curated HBV reference panel.
//      The bitscore-weighted majority vote of top hits assigns the genotype.
//
//      Confidence levels:
//        HIGH   – ≥80 % of total bitscore from one genotype
//        MEDIUM – 50–79 %
//        LOW    – <50 % or no hit above the % identity threshold
//
//      IMPORTANT: This is a BLAST pre-screen only.  For publication-quality
//      genotyping, run IQ-TREE manually on the alignments/ outputs.
//
// [VI] Gán genotype HBV (A–I) cho trình tự đồng thuận từng mẫu dùng BLAST.
//      Trình tự đồng thuận được BLAST với bộ tham chiếu HBV tuyển chọn.
//      Biểu quyết đa số theo trọng số bitscore của các hit hàng đầu gán genotype.
//
//      Mức độ tin cậy:
//        HIGH   – ≥80% tổng bitscore từ một genotype
//        MEDIUM – 50–79%
//        LOW    – <50% hoặc không có hit vượt ngưỡng % đồng nhất
//
//      QUAN TRỌNG: Đây chỉ là sàng lọc BLAST sơ bộ. Để gán genotype
//      chất lượng công bố, hãy chạy IQ-TREE thủ công trên đầu ra alignments/.
// ==============================================================================

process ASSIGN_GENOTYPE {

    tag "${sample_id}"
    label 'process_medium'

    // [EN] Save BLAST results and genotype TSV for audit
    // [VI] Lưu kết quả BLAST và TSV genotype để kiểm định
    publishDir "${params.outdir}/genotype/${sample_id}", mode: 'copy'

    input:
    // [EN] Consensus FASTA from BUILD_CONSENSUS
    // [VI] FASTA đồng thuận từ BUILD_CONSENSUS
    tuple val(sample_id), path(consensus_fasta)

    // [EN] HBV genotype reference FASTA panel
    // [VI] Bộ FASTA tham chiếu genotype HBV
    path reference_fasta

    output:
    // [EN] Genotype assignment TSV with confidence score
    // [VI] TSV gán genotype kèm mức độ tin cậy
    tuple val(sample_id), path("${sample_id}_genotype.tsv"),     emit: genotype_tsv

    // [EN] Raw BLAST tabular output (for manual inspection)
    // [VI] Đầu ra BLAST dạng bảng thô (để kiểm tra thủ công)
    path "${sample_id}_blast_raw.txt",                            emit: blast_raw

    script:
    """
    echo "[ASSIGN_GENOTYPE] ${sample_id}: Bắt đầu gán genotype / Starting genotype assignment..." >&2

    # [EN] Check if BLAST+ is installed
    # [VI] Kiểm tra BLAST+ đã được cài đặt
    if ! command -v blastn >/dev/null 2>&1; then
        echo "[ASSIGN_GENOTYPE] CẢNH BÁO / WARN: blastn không tìm thấy / not found." >&2
        echo "[ASSIGN_GENOTYPE] Gán UNKNOWN / Assigning UNKNOWN genotype." >&2
        echo "[ASSIGN_GENOTYPE] Cài đặt / Install: conda install -c bioconda blast" >&2
        touch "${sample_id}_blast_raw.txt"
        genotype_assign.py \\
            --blast-result "${sample_id}_blast_raw.txt" \\
            --sample-id    "${sample_id}" \\
            --min-pct-id   "${params.min_genotype_pct}" \\
            --out-tsv      "${sample_id}_genotype.tsv"
        exit 0
    fi

    # [EN] Build BLAST database from reference FASTA (per-task, isolated)
    # [VI] Xây dựng cơ sở dữ liệu BLAST từ FASTA tham chiếu (theo task, độc lập)
    echo "[ASSIGN_GENOTYPE] ${sample_id}: Xây dựng BLAST DB / Building BLAST DB..." >&2
    # [EN] Redirect both stdout/stderr to log file, consistent with blastn check above
    # [VI] Chuyển hướng cả stdout/stderr vào log file, nhất quán với kiểm tra blastn
    makeblastdb -in "${reference_fasta}" -dbtype nucl -out blast_refdb \\
        -title "HBV_genotype_refs" >"makeblastdb.log" 2>&1

    # [EN] Run BLAST: consensus vs reference panel
    # [VI] Chạy BLAST: đồng thuận so với bộ tham chiếu
    echo "[ASSIGN_GENOTYPE] ${sample_id}: Chạy blastn / Running blastn..." >&2
    blastn \\
        -query   "${consensus_fasta}" \\
        -db      blast_refdb \\
        -evalue  "${params.blast_evalue}" \\
        -outfmt  6 \\
        -max_hsps 1 \\
        -out     "${sample_id}_blast_raw.txt"

    # [EN] Assign genotype from BLAST results
    # [VI] Gán genotype từ kết quả BLAST
    genotype_assign.py \\
        --blast-result "${sample_id}_blast_raw.txt" \\
        --sample-id    "${sample_id}" \\
        --min-pct-id   "${params.min_genotype_pct}" \\
        --out-tsv      "${sample_id}_genotype.tsv"

    echo "[ASSIGN_GENOTYPE] ${sample_id}: Hoàn thành / Done." >&2
    """
}
