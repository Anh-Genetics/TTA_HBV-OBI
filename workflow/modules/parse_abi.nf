// ==============================================================================
// Module: PARSE_ABI
// [EN] Parse a single Sanger ABI (.ab1) chromatogram file into:
//      - FASTA file with the base-called sequence
//      - TSV file with per-position Phred quality scores and IUPAC flags
//
// [VI] Phân tích một tệp chromatogram Sanger ABI (.ab1) thành:
//      - Tệp FASTA chứa trình tự base-called
//      - Tệp TSV ghi điểm chất lượng Phred từng vị trí và cờ IUPAC
//
// [EN] Uses Biopython for reliable ABI channel parsing.
// [VI] Sử dụng Biopython để phân tích kênh ABI một cách tin cậy.
// ==============================================================================

process PARSE_ABI {

    tag "${sample_id}:${direction}"
    label 'process_low'

    // [EN] Save raw read FASTA and quality TSV for inspection
    // [VI] Lưu FASTA đoạn đọc thô và TSV chất lượng để kiểm tra
    publishDir "${params.outdir}/reads/raw/${sample_id}", mode: 'copy'

    input:
    // [EN] Tuple: (sample_id, direction, ab1_file)
    // [VI] Bộ ba: (mã_mẫu, hướng_đọc, tệp_ab1)
    tuple val(sample_id), val(direction), path(ab1_file)

    output:
    // [EN] FASTA + quality TSV, keyed by sample_id + direction
    // [VI] FASTA và TSV chất lượng, theo mã mẫu và hướng đọc
    tuple val(sample_id), val(direction), path("${sample_id}_${direction}.fasta"),          emit: read_fasta
    tuple val(sample_id), val(direction), path("${sample_id}_${direction}_quality.tsv"),    emit: quality_tsv

    script:
    """
    echo "[PARSE_ABI] ${sample_id} (${direction}): Phân tích tệp ABI / Parsing ABI file: ${ab1_file}" >&2

    parse_abi.py \\
        --input      "${ab1_file}" \\
        --sample-id  "${sample_id}" \\
        --direction  "${direction}" \\
        --out-fasta  "${sample_id}_${direction}.fasta" \\
        --out-qual   "${sample_id}_${direction}_quality.tsv"

    echo "[PARSE_ABI] ${sample_id} (${direction}): Hoàn thành / Done." >&2
    """
}
