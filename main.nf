#!/usr/bin/env nextflow
// ==============================================================================
// TTA_HBV-OBI  –  HBV pre-S1/pre-S2/S Sanger Analysis Pipeline
// Nextflow DSL2 | Phiên bản / Version: 0.1.0-MVP
//
// [EN] PURPOSE
//   Analyse Sanger sequencing reads (ABI .ab1 files) covering the HBV
//   pre-S1 / pre-S2 / S region to:
//     1. Parse and quality-trim ABI chromatograms
//     2. Build per-sample consensus sequences
//     3. Assign HBV genotype (A-I) via BLAST pre-screen
//     4. Annotate variants against OBI-associated mutation catalogue
//     5. Generate structured TSV + HTML summary reports
//
// [VI] MỤC ĐÍCH
//   Phân tích các đoạn đọc Sanger (tệp ABI .ab1) bao phủ vùng
//   pre-S1 / pre-S2 / S của HBV để:
//     1. Phân tích và cắt tỉa chất lượng chromatogram ABI
//     2. Xây dựng trình tự đồng thuận cho từng mẫu
//     3. Gán genotype HBV (A-I) qua sàng lọc BLAST sơ bộ
//     4. Chú thích biến thể theo danh mục đột biến liên quan OBI
//     5. Tạo báo cáo tóm tắt dạng TSV + HTML
//
// [EN] IMPORTANT LIMITATIONS (MVP):
//   - Genotype assignment uses BLAST pre-screen only; full ML tree is scaffolded
//     but not run automatically.
//   - OBI variant interpretation is based on published literature catalogue;
//     it is NOT a clinical diagnostic tool.
//   - Minor variants / quasispecies within a single Sanger read cannot be
//     reliably detected; dual-peak positions are flagged as IUPAC ambiguities.
//
// [VI] HẠN CHẾ QUAN TRỌNG (MVP):
//   - Gán genotype chỉ dùng sàng lọc BLAST; cây ML không được dựng tự động.
//   - Diễn giải biến thể OBI dựa trên danh mục tài liệu đã công bố;
//     KHÔNG phải công cụ chẩn đoán lâm sàng.
//   - Không thể phát hiện minor variant / quasispecies từ một đoạn đọc Sanger
//     đơn lẻ; các vị trí peak đôi được ghi nhận là mã IUPAC mơ hồ.
//
// [EN] USAGE:
//   nextflow run main.nf --sample_sheet data/example/sample_sheet.csv \
//                        --outdir results/ \
//                        -profile local
//
// [VI] CÁCH DÙNG:
//   nextflow run main.nf --sample_sheet data/example/sample_sheet.csv \
//                        --outdir results/ \
//                        -profile local
// ==============================================================================

nextflow.enable.dsl = 2

// ─── Tham số pipeline / Pipeline parameters ──────────────────────────────────
params.sample_sheet     = null
params.outdir           = "results"
params.reference_fasta  = "${projectDir}/data/references/hbv_genotype_refs.fasta"
params.mutation_db      = "${projectDir}/data/references/obi_mutation_catalogue.tsv"

// Cắt tỉa chất lượng / QC & trimming
params.min_quality      = 20     // ngưỡng chất lượng Phred / Phred quality threshold
params.min_length       = 200    // độ dài tối thiểu sau cắt (bp) / min length after trim (bp)
params.trim_ends_bases  = 20     // số base cố định cắt mỗi đầu / fixed bases to trim each end

// Xây dựng đồng thuận / Consensus building
params.min_overlap      = 80     // vùng chồng tối thiểu (bp) / minimum overlap (bp)
params.conflict_policy  = "iupac"  // chính sách xử lý xung đột / conflict policy: iupac|majority|n

// Gán genotype / Genotype assignment
params.blast_evalue     = "1e-10"
params.min_genotype_pct = 90.0   // % đồng nhất tối thiểu / minimum % identity

// Báo cáo / Report
params.report_format    = "html"   // định dạng / format: html|tsv|markdown

// ─── Nạp module / Include modules ────────────────────────────────────────────
include { VALIDATE_SAMPLESHEET  } from './workflow/modules/validate_inputs.nf'
include { PARSE_ABI             } from './workflow/modules/parse_abi.nf'
include { TRIM_READS            } from './workflow/modules/trim_reads.nf'
include { BUILD_CONSENSUS       } from './workflow/modules/build_consensus.nf'
include { ALIGN_CONSENSUS       } from './workflow/modules/align_consensus.nf'
include { ASSIGN_GENOTYPE       } from './workflow/modules/genotype_assign.nf'
include { ANNOTATE_VARIANTS     } from './workflow/modules/annotate_variants.nf'
include { GENERATE_REPORT       } from './workflow/modules/generate_report.nf'

// ─── Workflow chính / Main workflow ──────────────────────────────────────────
workflow {

    // ── Bước 0: Kiểm tra tham số bắt buộc / Step 0: Check required params ──
    if (!params.sample_sheet) {
        error """
        LỖI / ERROR: --sample_sheet là bắt buộc / is required.

        [VI] Cách dùng:
          nextflow run main.nf \\
            --sample_sheet data/example/sample_sheet.csv \\
            --outdir results/

        [EN] Usage:
          nextflow run main.nf \\
            --sample_sheet data/example/sample_sheet.csv \\
            --outdir results/

        Xem / See README.md để biết thêm chi tiết / for full usage instructions.
        """.stripIndent()
    }

    // In thông tin khởi động / Print startup banner
    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║       TTA_HBV-OBI  –  Sanger Analysis Pipeline              ║
    ║       Phiên bản / Version: 0.1.0-MVP                        ║
    ╠══════════════════════════════════════════════════════════════╣
    ║  [VI] Pipeline phân tích Sanger HBV pre-S1/pre-S2/S         ║
    ║  [EN] HBV pre-S1/pre-S2/S Sanger sequencing pipeline        ║
    ╚══════════════════════════════════════════════════════════════╝

    Tệp mô tả mẫu / Sample sheet : ${params.sample_sheet}
    Thư mục đầu ra / Output dir  : ${params.outdir}
    Tham chiếu / Reference       : ${params.reference_fasta}
    Danh mục đột biến / Mut. DB  : ${params.mutation_db}
    Chính sách xung đột / Conflict: ${params.conflict_policy}
    Định dạng báo cáo / Report   : ${params.report_format}
    """.stripIndent()

    // ── Bước 1: Kiểm tra tệp mô tả mẫu / Step 1: Validate sample sheet ─────
    log.info "[TTA_HBV-OBI] Bước 1/8: Kiểm tra tệp mô tả mẫu / Step 1/8: Validating sample sheet..."
    sample_sheet_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
    VALIDATE_SAMPLESHEET(sample_sheet_ch)

    // Phân tích TSV đã kiểm tra / Parse validated TSV
    // Định dạng: sample_id \t ab1_forward \t ab1_reverse (hoặc NONE)
    validated_rows_ch = VALIDATE_SAMPLESHEET.out.validated_tsv
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def fwd = file(row.ab1_forward, checkIfExists: true)
            def rev = (row.ab1_reverse == 'NONE') ? null : file(row.ab1_reverse, checkIfExists: true)
            tuple(row.sample_id, fwd, rev)
        }

    // ── Bước 2: Phân tích tệp ABI / Step 2: Parse ABI files ─────────────────
    log.info "[TTA_HBV-OBI] Bước 2/8: Phân tích tệp ABI / Step 2/8: Parsing ABI files..."

    // Tách kênh đọc xuôi và ngược / Split forward and reverse read channels
    fwd_reads_ch = validated_rows_ch.map { sid, fwd, rev -> tuple(sid, 'forward', fwd) }
    rev_reads_ch = validated_rows_ch
        .filter { sid, fwd, rev -> rev != null }
        .map     { sid, fwd, rev -> tuple(sid, 'reverse', rev) }

    all_reads_ch = fwd_reads_ch.mix(rev_reads_ch)
    PARSE_ABI(all_reads_ch)

    // ── Bước 3: Cắt tỉa chất lượng / Step 3: Quality trimming ───────────────
    log.info "[TTA_HBV-OBI] Bước 3/8: Cắt tỉa chất lượng / Step 3/8: Quality trimming reads..."
    TRIM_READS(PARSE_ABI.out.read_fasta, PARSE_ABI.out.quality_tsv)

    // ── Bước 4: Nhóm và xây dựng đồng thuận / Step 4: Group & build consensus
    // [EN] groupTuple() ensures ALL reads for a sample are collected before
    //      consensus building – prevents race conditions.
    // [VI] groupTuple() đảm bảo TẤT CẢ đoạn đọc của một mẫu được thu thập
    //      trước khi xây dựng đồng thuận – tránh xung đột tiến trình song song.
    log.info "[TTA_HBV-OBI] Bước 4/8: Xây dựng trình tự đồng thuận / Step 4/8: Building consensus sequences..."
    grouped_reads_ch = TRIM_READS.out.trimmed_fasta
        .map    { sid, dir, fasta, stats -> tuple(sid, fasta) }
        .groupTuple(by: 0)

    BUILD_CONSENSUS(grouped_reads_ch)

    // ── Bước 5: Căn chỉnh với tham chiếu / Step 5: Align to references ──────
    log.info "[TTA_HBV-OBI] Bước 5/8: Căn chỉnh với tham chiếu HBV / Step 5/8: Aligning to HBV references (MAFFT)..."
    ref_ch = Channel.fromPath(params.reference_fasta)
    ALIGN_CONSENSUS(BUILD_CONSENSUS.out.consensus_fasta, ref_ch)

    // ── Bước 6: Gán genotype / Step 6: Assign genotype ──────────────────────
    log.info "[TTA_HBV-OBI] Bước 6/8: Gán genotype HBV / Step 6/8: Assigning HBV genotypes (BLAST)..."
    ASSIGN_GENOTYPE(BUILD_CONSENSUS.out.consensus_fasta, ref_ch)

    // ── Bước 7: Chú thích biến thể / Step 7: Annotate variants ──────────────
    log.info "[TTA_HBV-OBI] Bước 7/8: Chú thích biến thể OBI / Step 7/8: Annotating OBI variants..."
    mut_db_ch = Channel.fromPath(params.mutation_db)
    ANNOTATE_VARIANTS(
        BUILD_CONSENSUS.out.consensus_fasta,
        ALIGN_CONSENSUS.out.alignment_fasta,
        mut_db_ch
    )

    // ── Bước 8: Tạo báo cáo / Step 8: Generate reports ──────────────────────
    // [EN] Collect all per-sample outputs by joining on sample_id before report
    //      generation – no shared file writes during parallel execution.
    // [VI] Thu thập tất cả đầu ra theo mẫu bằng join trước khi tạo báo cáo
    //      – không ghi tệp dùng chung trong quá trình thực thi song song.
    log.info "[TTA_HBV-OBI] Bước 8/8: Tạo báo cáo / Step 8/8: Generating per-sample reports..."

    report_inputs_ch = BUILD_CONSENSUS.out.consensus_fasta
        .join(ASSIGN_GENOTYPE.out.genotype_tsv,     by: 0)
        .join(ANNOTATE_VARIANTS.out.annotation_tsv, by: 0)
        .join(TRIM_READS.out.trim_summary,          by: 0)
        .join(BUILD_CONSENSUS.out.consensus_stats,  by: 0)

    GENERATE_REPORT(report_inputs_ch)

    // Chờ tất cả báo cáo hoàn thành và in tóm tắt
    // Wait for all reports to finish and print summary
    GENERATE_REPORT.out.sample_report.collect().view { reports ->
        log.info "\n[TTA_HBV-OBI] [OK] Pipeline hoàn thành / Pipeline complete."
        log.info "[TTA_HBV-OBI]    ${reports.size()} mẫu / sample(s) đã xử lý / processed."
        log.info "[TTA_HBV-OBI]    Kết quả lưu tại / Results in: ${params.outdir}"
    }
}

// ─── Sự kiện hoàn thành pipeline / Pipeline completion event ────────────────
workflow.onComplete {
    def status = workflow.success ? "THANH CONG / SUCCESS [OK]" : "THAT BAI / FAILED [FAIL]"
    log.info """
    ──────────────────────────────────────────────────────────────
    [TTA_HBV-OBI] Kết quả / Result  : ${status}
    [TTA_HBV-OBI] Thời gian / Duration: ${workflow.duration}
    [TTA_HBV-OBI] Mã thoát / Exit code: ${workflow.exitStatus}
    [TTA_HBV-OBI] Kết quả tại / Output: ${params.outdir}
    ──────────────────────────────────────────────────────────────
    """.stripIndent()
}

// ─── Sự kiện lỗi pipeline / Pipeline error event ────────────────────────────
workflow.onError {
    log.error "[TTA_HBV-OBI] LỖI PIPELINE / Pipeline error: ${workflow.errorMessage}"
}
