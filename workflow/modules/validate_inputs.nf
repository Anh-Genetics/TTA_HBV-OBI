// ==============================================================================
// Module: VALIDATE_SAMPLESHEET
// [EN] Validate the CSV sample sheet before any analysis begins.
//      Checks: required columns present, file paths exist, no duplicate sample_id.
//      Writes a validated TSV to avoid repeated validation downstream.
// [VI] Kiểm tra tệp CSV mô tả mẫu trước khi bất kỳ bước phân tích nào bắt đầu.
//      Kiểm tra: cột bắt buộc có mặt, đường dẫn tệp tồn tại, không trùng mã mẫu.
//      Ghi TSV đã kiểm tra để tránh kiểm tra lại ở các bước sau.
// ==============================================================================

process VALIDATE_SAMPLESHEET {

    tag "validate_samplesheet"
    label 'process_low'

    // [EN] Publish the validation report alongside the validated TSV
    // [VI] Lưu báo cáo kiểm tra cùng với TSV đã xác thực
    publishDir "${params.outdir}/qc/samplesheet", mode: 'copy'

    input:
    // [EN] Raw sample sheet CSV from the user
    // [VI] Tệp CSV mô tả mẫu gốc từ người dùng
    path sample_sheet

    output:
    // [EN] Validated TSV used by downstream processes
    // [VI] TSV đã kiểm tra dùng cho các bước xử lý tiếp theo
    path "validated_samplesheet.tsv", emit: validated_tsv

    // [EN] Human-readable validation report for audit purposes
    // [VI] Báo cáo kiểm tra dạng văn bản để đánh giá và kiểm định
    path "samplesheet_validation_report.txt", emit: validation_report

    script:
    """
    # [EN] Informational: show which sample sheet is being validated
    # [VI] Thông báo: hiển thị tệp mô tả mẫu đang được kiểm tra
    echo "[VALIDATE_SAMPLESHEET] Bắt đầu kiểm tra / Starting validation: ${sample_sheet}" >&2

    validate_samplesheet.py \\
        --input   "${sample_sheet}" \\
        --output  validated_samplesheet.tsv \\
        --report  samplesheet_validation_report.txt

    # [EN] Report validation outcome to the console
    # [VI] Báo cáo kết quả kiểm tra ra màn hình
    echo "[VALIDATE_SAMPLESHEET] Hoàn thành / Complete. Xem / See: samplesheet_validation_report.txt" >&2
    """
}
