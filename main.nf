#!/usr/bin/env nextflow
// ==============================================================================
// TTA_HBV-OBI  –  HBV pre-S1/pre-S2/S Sanger Analysis Pipeline
// Nextflow DSL2 | Version: 0.1.0-MVP
// ==============================================================================
//
// PURPOSE
//   Analyse Sanger sequencing reads (ABI .ab1 files) covering the HBV
//   pre-S1 / pre-S2 / S region to:
//     1. Parse and quality-trim ABI chromatograms
//     2. Build per-sample consensus sequences
//     3. Assign HBV genotype (A–I) via alignment + phylogenetic pre-screen
//     4. Annotate variants against OBI-associated mutation catalogue
//     5. Generate a structured TSV + HTML summary report
//
// IMPORTANT LIMITATIONS (MVP)
//   - Genotype assignment uses alignment-based pre-screen only; a full
//     maximum-likelihood tree is scaffolded but not run automatically to
//     avoid long runtimes on a laptop.
//   - OBI variant interpretation is based on a curated literature catalogue;
//     it is NOT a clinical diagnosis tool.
//   - Minor variants / quasispecies within a single Sanger read cannot be
//     reliably detected; dual-peak positions are flagged as IUPAC ambiguities.
//
// USAGE
//   nextflow run main.nf --sample_sheet data/example/sample_sheet.csv \
//                        --outdir results/ \
//                        -profile local
//
// ==============================================================================

nextflow.enable.dsl = 2

// --------------------------------------------------------------------------- //
// Pipeline parameters (overridable on CLI or in nextflow.config)
// --------------------------------------------------------------------------- //
params.sample_sheet     = null
params.outdir           = "results"
params.reference_fasta  = "${projectDir}/data/references/hbv_genotype_refs.fasta"
params.mutation_db      = "${projectDir}/data/references/obi_mutation_catalogue.tsv"

// QC / trimming
params.min_quality      = 20     // phred-like quality threshold for base trimming
params.min_length       = 200    // minimum sequence length after trimming (bp)
params.trim_ends_bases  = 20     // always trim this many bases from each end

// Consensus building
params.min_overlap      = 80     // minimum overlap (bp) required to merge F/R reads
params.conflict_policy  = "iupac"  // "iupac" | "majority" | "n"

// Genotype assignment
params.blast_evalue     = "1e-10"
params.min_genotype_pct = 90.0   // minimum % identity to call a genotype

// Report
params.report_format    = "html"   // "html" | "tsv" | "markdown"

// --------------------------------------------------------------------------- //
// Include modules
// --------------------------------------------------------------------------- //
include { VALIDATE_SAMPLESHEET  } from './workflow/modules/validate_inputs.nf'
include { PARSE_ABI             } from './workflow/modules/parse_abi.nf'
include { TRIM_READS            } from './workflow/modules/trim_reads.nf'
include { BUILD_CONSENSUS       } from './workflow/modules/build_consensus.nf'
include { ALIGN_CONSENSUS       } from './workflow/modules/align_consensus.nf'
include { ASSIGN_GENOTYPE       } from './workflow/modules/genotype_assign.nf'
include { ANNOTATE_VARIANTS     } from './workflow/modules/annotate_variants.nf'
include { GENERATE_REPORT       } from './workflow/modules/generate_report.nf'

// --------------------------------------------------------------------------- //
// Helper: emit a warning to log
// --------------------------------------------------------------------------- //
def warnLog(String msg) {
    log.warn "  [TTA_HBV-OBI] ${msg}"
}

// --------------------------------------------------------------------------- //
// Main workflow
// --------------------------------------------------------------------------- //
workflow {

    // ---------------------------------------------------------------------- //
    // 0. Startup checks
    // ---------------------------------------------------------------------- //
    if (!params.sample_sheet) {
        error """
        ERROR: --sample_sheet is required.
        Usage:
          nextflow run main.nf \\
            --sample_sheet data/example/sample_sheet.csv \\
            --outdir results/

        See README.md for full usage instructions.
        """.stripIndent()
    }

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║         TTA_HBV-OBI  Sanger Analysis Pipeline        ║
    ║         Version: 0.1.0-MVP                           ║
    ╚══════════════════════════════════════════════════════╝
    Sample sheet  : ${params.sample_sheet}
    Output dir    : ${params.outdir}
    Reference     : ${params.reference_fasta}
    Mutation DB   : ${params.mutation_db}
    Conflict      : ${params.conflict_policy}
    Report format : ${params.report_format}
    """.stripIndent()

    // ---------------------------------------------------------------------- //
    // 1. Validate sample sheet – emits channel: [sample_id, fwd_ab1, rev_ab1]
    //    rev_ab1 may be the string "NONE" for single-read samples
    // ---------------------------------------------------------------------- //
    sample_sheet_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)

    VALIDATE_SAMPLESHEET(sample_sheet_ch)

    // Parse the validated TSV emitted by the process (one row per sample)
    // Format: sample_id \t ab1_forward \t ab1_reverse (or NONE)
    validated_rows_ch = VALIDATE_SAMPLESHEET.out.validated_tsv
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def fwd = file(row.ab1_forward, checkIfExists: true)
            def rev = (row.ab1_reverse == 'NONE') ? null : file(row.ab1_reverse, checkIfExists: true)
            tuple(row.sample_id, fwd, rev)
        }

    // ---------------------------------------------------------------------- //
    // 2. Parse ABI files  →  per-read FASTA + quality TSV
    // ---------------------------------------------------------------------- //
    // Emit one tuple per read; reads without a reverse get a single entry
    fwd_reads_ch = validated_rows_ch.map { sid, fwd, rev -> tuple(sid, 'forward', fwd) }
    rev_reads_ch = validated_rows_ch
        .filter { sid, fwd, rev -> rev != null }
        .map     { sid, fwd, rev -> tuple(sid, 'reverse', rev) }

    all_reads_ch = fwd_reads_ch.mix(rev_reads_ch)

    PARSE_ABI(all_reads_ch)

    // ---------------------------------------------------------------------- //
    // 3. Trim reads by quality
    // ---------------------------------------------------------------------- //
    TRIM_READS(PARSE_ABI.out.read_fasta, PARSE_ABI.out.quality_tsv)

    // ---------------------------------------------------------------------- //
    // 4. Group trimmed reads by sample_id, then build consensus
    //    This grouping step avoids any race conditions: each sample's reads
    //    are collected before the consensus step begins.
    // ---------------------------------------------------------------------- //
    // TRIM_READS emits: [sample_id, direction, trimmed_fasta, trim_stats]
    grouped_reads_ch = TRIM_READS.out.trimmed_fasta
        .map    { sid, dir, fasta, stats -> tuple(sid, fasta) }
        .groupTuple(by: 0)   // group all reads for a given sample_id together

    BUILD_CONSENSUS(grouped_reads_ch)

    // ---------------------------------------------------------------------- //
    // 5. Align consensus to reference (MAFFT)
    // ---------------------------------------------------------------------- //
    ref_ch = Channel.fromPath(params.reference_fasta)

    ALIGN_CONSENSUS(BUILD_CONSENSUS.out.consensus_fasta, ref_ch)

    // ---------------------------------------------------------------------- //
    // 6. Assign genotype
    // ---------------------------------------------------------------------- //
    ASSIGN_GENOTYPE(BUILD_CONSENSUS.out.consensus_fasta, ref_ch)

    // ---------------------------------------------------------------------- //
    // 7. Annotate variants against OBI mutation catalogue
    // ---------------------------------------------------------------------- //
    mut_db_ch = Channel.fromPath(params.mutation_db)

    ANNOTATE_VARIANTS(
        BUILD_CONSENSUS.out.consensus_fasta,
        ALIGN_CONSENSUS.out.alignment_fasta,
        mut_db_ch
    )

    // ---------------------------------------------------------------------- //
    // 8. Generate per-sample + run-level reports
    //    Collect all per-sample outputs before writing the summary to avoid
    //    partial/race-condition writes.
    // ---------------------------------------------------------------------- //
    // Combine per-sample outputs keyed on sample_id
    report_inputs_ch = BUILD_CONSENSUS.out.consensus_fasta
        .join(ASSIGN_GENOTYPE.out.genotype_tsv,     by: 0)
        .join(ANNOTATE_VARIANTS.out.annotation_tsv, by: 0)
        .join(TRIM_READS.out.trim_summary,          by: 0)
        .join(BUILD_CONSENSUS.out.consensus_stats,  by: 0)

    GENERATE_REPORT(report_inputs_ch)

    // Collect the run-level summary TSV by waiting for ALL per-sample reports
    all_reports_ch = GENERATE_REPORT.out.sample_report.collect()
    GENERATE_REPORT.out.sample_report.collect().view { reports ->
        log.info "\n[TTA_HBV-OBI] Pipeline complete. ${reports.size()} sample(s) processed."
        log.info "[TTA_HBV-OBI] Results written to: ${params.outdir}"
    }
}

// --------------------------------------------------------------------------- //
// Workflow completion handler
// --------------------------------------------------------------------------- //
workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    log.info """
    ──────────────────────────────────────────────
    Pipeline ${status}
    Duration : ${workflow.duration}
    Exit code: ${workflow.exitStatus}
    Results  : ${params.outdir}
    ──────────────────────────────────────────────
    """.stripIndent()
}

workflow.onError {
    log.error "[TTA_HBV-OBI] Pipeline error: ${workflow.errorMessage}"
}
