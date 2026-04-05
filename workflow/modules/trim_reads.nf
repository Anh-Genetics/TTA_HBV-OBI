// workflow/modules/trim_reads.nf
// Quality-trims a parsed Sanger read FASTA using the associated quality TSV.
// Inputs:  read_fasta  [sample_id, direction, fasta_file]
//          quality_tsv [sample_id, direction, tsv_file]
// Outputs: trimmed_fasta [sample_id, direction, trimmed_fasta, stats_tsv]
//          trim_summary  [sample_id, stats_tsv]   (joined by sample_id for report)

process TRIM_READS {
    label 'low'
    tag   "${sample_id}_${direction}"

    publishDir "${params.outdir}/trimmed_reads/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(direction), path(fasta_file)
    tuple val(sample_id), val(direction), path(quality_tsv)

    output:
    tuple val(sample_id), val(direction), path("${sample_id}_${direction}_trimmed.fasta"), path("${sample_id}_${direction}_trim_stats.tsv"), emit: trimmed_fasta
    tuple val(sample_id), path("${sample_id}_${direction}_trim_stats.tsv"), emit: trim_summary

    script:
    """
    trim_reads.py \\
        --fasta         "${fasta_file}" \\
        --quality       "${quality_tsv}" \\
        --sample-id     "${sample_id}" \\
        --direction     "${direction}" \\
        --min-quality   ${params.min_quality} \\
        --min-length    ${params.min_length} \\
        --trim-ends     ${params.trim_ends_bases} \\
        --out-fasta     "${sample_id}_${direction}_trimmed.fasta" \\
        --out-stats     "${sample_id}_${direction}_trim_stats.tsv"
    """
}
