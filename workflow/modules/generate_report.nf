// workflow/modules/generate_report.nf
// Generates a per-sample report (HTML, TSV, or Markdown) summarising:
//   - QC / trim statistics
//   - Consensus quality
//   - Genotype assignment
//   - Variant annotation with OBI relevance classification
//
// Each sample runs in its own isolated process → no shared-file race conditions.
// The run-level summary is assembled from the collect()ed per-sample outputs
// in main.nf (outside this process).
//
// Inputs:  [sample_id, consensus_fasta, genotype_tsv, annotation_tsv, trim_summary, consensus_stats]
// Outputs: sample_report [sample_id, report_file]

process GENERATE_REPORT {
    label 'low'
    tag   "${sample_id}"

    publishDir "${params.outdir}/reports/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id),
          path(consensus_fasta),
          path(genotype_tsv),
          path(annotation_tsv),
          path(trim_summary),
          path(consensus_stats)

    output:
    tuple val(sample_id), path("${sample_id}_report.${params.report_format}"), emit: sample_report

    script:
    """
    generate_report.py \\
        --sample-id     "${sample_id}" \\
        --consensus     "${consensus_fasta}" \\
        --genotype      "${genotype_tsv}" \\
        --variants      "${annotation_tsv}" \\
        --trim-stats    "${trim_summary}" \\
        --consensus-stats "${consensus_stats}" \\
        --format        "${params.report_format}" \\
        --out           "${sample_id}_report.${params.report_format}"
    """
}
