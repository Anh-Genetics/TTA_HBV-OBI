// workflow/modules/parse_abi.nf
// Parses a single ABI .ab1 file into a FASTA sequence and quality TSV.
// Inputs: [sample_id, direction, ab1_file]
// Outputs: [sample_id, direction, read_fasta], [sample_id, direction, quality_tsv]

process PARSE_ABI {
    label 'low'
    tag   "${sample_id}_${direction}"

    publishDir "${params.outdir}/parsed_reads/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(direction), path(ab1_file)

    output:
    tuple val(sample_id), val(direction), path("${sample_id}_${direction}.fasta"), emit: read_fasta
    tuple val(sample_id), val(direction), path("${sample_id}_${direction}_quality.tsv"), emit: quality_tsv

    script:
    """
    parse_abi.py \\
        --input     "${ab1_file}" \\
        --sample-id "${sample_id}" \\
        --direction "${direction}" \\
        --out-fasta "${sample_id}_${direction}.fasta" \\
        --out-qual  "${sample_id}_${direction}_quality.tsv"
    """
}
