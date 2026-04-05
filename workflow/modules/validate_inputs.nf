// workflow/modules/validate_inputs.nf
// Validates the sample sheet CSV and confirms that all referenced files exist.
// Emits a clean TSV for downstream parsing.

process VALIDATE_SAMPLESHEET {
    label 'low'
    tag   "validate_samplesheet"

    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path sample_sheet

    output:
    path "validated_samplesheet.tsv",  emit: validated_tsv
    path "validation_report.txt",      emit: validation_report

    script:
    """
    validate_samplesheet.py \\
        --input   "${sample_sheet}" \\
        --output  validated_samplesheet.tsv \\
        --report  validation_report.txt
    """
}
