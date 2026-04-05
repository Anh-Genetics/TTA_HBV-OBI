"""
Tests for the TTA_HBV-OBI pipeline Python helper scripts.

These tests focus on:
  - validate_samplesheet.py: CSV parsing, file validation logic
  - build_consensus.py: overlap detection, IUPAC merging, RC
  - genotype_assign.py: BLAST result parsing, genotype voting
  - annotate_variants.py: variant detection, catalogue lookup
  - trim_reads.py: quality trimming logic
  - generate_report.py: report generation

All tests use synthetic in-memory data; no real .ab1 files required.
"""

import csv
import io
import sys
import os
import tempfile
import textwrap
from pathlib import Path

import pytest

# Add workflow/bin to path so we can import the scripts
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "bin"))


# ============================================================================
# Helper utilities
# ============================================================================

def write_temp_file(content: str, suffix: str = ".txt") -> Path:
    """Write content to a temporary file and return the path."""
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=suffix, delete=False, encoding="utf-8"
    )
    tmp.write(content)
    tmp.close()
    return Path(tmp.name)


def write_temp_fasta(header: str, seq: str) -> Path:
    return write_temp_file(f">{header}\n{seq}\n", suffix=".fasta")


def write_temp_tsv(rows: list[dict], fieldnames: list[str]) -> Path:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return write_temp_file(buf.getvalue(), suffix=".tsv")


# ============================================================================
# Tests: build_consensus module
# ============================================================================

class TestBuildConsensus:
    """Tests for consensus-building logic."""

    def _import_bc(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "build_consensus",
            Path(__file__).parent.parent / "workflow" / "bin" / "build_consensus.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_reverse_complement_basic(self):
        bc = self._import_bc()
        assert bc.reverse_complement("ATCG") == "CGAT"
        assert bc.reverse_complement("AAAA") == "TTTT"
        assert bc.reverse_complement("GCGC") == "GCGC"

    def test_reverse_complement_iupac(self):
        bc = self._import_bc()
        # R (A/G) → Y (C/T)
        rc = bc.reverse_complement("R")
        assert rc == "Y"

    def test_iupac_merge_same_base(self):
        bc = self._import_bc()
        assert bc.iupac_merge("A", "A") == "A"
        assert bc.iupac_merge("C", "C") == "C"

    def test_iupac_merge_ag(self):
        bc = self._import_bc()
        assert bc.iupac_merge("A", "G") == "R"
        assert bc.iupac_merge("G", "A") == "R"

    def test_iupac_merge_unknown_pair(self):
        bc = self._import_bc()
        # A + C + G not a simple pair – should return N
        assert bc.iupac_merge("A", "N") == "N"

    def test_find_overlap_exact(self):
        bc = self._import_bc()
        fwd = "ATCGATCGATCG" + "GCTAGCTAGCTA"
        rev_rc = "GCTAGCTAGCTA" + "TTTTTTTTTTTT"
        ol = bc.find_overlap(fwd, rev_rc, min_overlap=10)
        assert ol == 12  # exact 12-bp overlap

    def test_find_overlap_no_overlap(self):
        bc = self._import_bc()
        fwd = "AAAAAAAAAAAAAAAAAAAAAA"
        rev_rc = "CCCCCCCCCCCCCCCCCCCCC"
        ol = bc.find_overlap(fwd, rev_rc, min_overlap=10)
        assert ol == 0

    def test_merge_overlap_iupac(self):
        bc = self._import_bc()
        # fwd last 6 bp = AATTGG
        # rev_rc first 6 bp = AACTGG  (conflict at pos 3: T vs C → Y)
        result = bc.merge_overlap("XXXXXAATTGG", "AACTGGYYYYY", overlap=6, conflict="iupac")
        # unique_fwd: XXXXX, overlap merged, unique_rev: YYYYY
        assert result.startswith("XXXXX")
        assert result.endswith("YYYYY")
        overlap_part = result[5:11]
        assert overlap_part[2] == "Y"  # T+C → Y

    def test_merge_overlap_n_policy(self):
        bc = self._import_bc()
        result = bc.merge_overlap("XXXXXAATTGG", "AACTGGYYYYY", overlap=6, conflict="n")
        overlap_part = result[5:11]
        assert overlap_part[2] == "N"

    def test_is_pass_read(self):
        bc = self._import_bc()
        assert bc.is_pass_read("sample_forward_trimmed status=PASS", "ATCG") is True
        assert bc.is_pass_read("sample_forward_trimmed status=FAIL", "ATCG") is False
        assert bc.is_pass_read("header", "") is False

    def test_full_single_read_consensus(self, tmp_path):
        bc = self._import_bc()
        fasta = tmp_path / "fwd_trimmed.fasta"
        fasta.write_text(">SAMP001_forward_trimmed status=PASS\nATCGATCGATCG\n")

        out_fasta = str(tmp_path / "consensus.fasta")
        out_stats = str(tmp_path / "consensus_stats.tsv")

        # Patch sys.argv and call main
        old_argv = sys.argv
        sys.argv = [
            "build_consensus.py",
            "--sample-id", "SAMP001",
            "--reads", str(fasta),
            "--min-overlap", "10",
            "--conflict", "iupac",
            "--out-fasta", out_fasta,
            "--out-stats", out_stats,
        ]
        try:
            bc.main()
        finally:
            sys.argv = old_argv

        result = Path(out_fasta).read_text()
        assert ">SAMP001" in result
        assert "ATCGATCGATCG" in result


# ============================================================================
# Tests: validate_samplesheet module
# ============================================================================

class TestValidateSamplesheet:
    """Tests for sample sheet validation logic."""

    def _import_vs(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "validate_samplesheet",
            Path(__file__).parent.parent / "workflow" / "bin" / "validate_samplesheet.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_normalise_header(self):
        vs = self._import_vs()
        result = vs.normalise_header(["Sample_ID", " Ab1_Forward "])
        assert "sample_id" in result
        assert "ab1_forward" in result

    def test_validate_file_path_none(self):
        vs = self._import_vs()
        path, errors = vs.validate_file_path("", "ab1_forward", "SAMP001")
        assert path == "NONE"
        assert errors == []

    def test_validate_file_path_missing(self, tmp_path):
        vs = self._import_vs()
        missing = str(tmp_path / "missing.ab1")
        path, errors = vs.validate_file_path(missing, "ab1_forward", "SAMP001")
        assert any("[ERROR]" in e for e in errors)

    def test_validate_file_path_exists(self, tmp_path):
        vs = self._import_vs()
        ab1 = tmp_path / "sample.ab1"
        ab1.write_bytes(b"")  # empty but exists
        path, errors = vs.validate_file_path(str(ab1), "ab1_forward", "SAMP001")
        assert errors == []
        assert str(ab1) in path

    def test_validate_file_path_wrong_ext(self, tmp_path):
        vs = self._import_vs()
        txt = tmp_path / "sample.txt"
        txt.write_bytes(b"")
        path, errors = vs.validate_file_path(str(txt), "ab1_forward", "SAMP001")
        assert any("[WARN]" in e for e in errors)

    def test_full_validation_pass(self, tmp_path):
        vs = self._import_vs()
        # Create dummy ab1 files
        fwd = tmp_path / "S001_F.ab1"
        rev = tmp_path / "S001_R.ab1"
        fwd.write_bytes(b"")
        rev.write_bytes(b"")

        csv_content = f"sample_id,ab1_forward,ab1_reverse\nS001,{fwd},{rev}\n"
        sheet = tmp_path / "sheet.csv"
        sheet.write_text(csv_content)
        out_tsv = str(tmp_path / "validated.tsv")
        report  = str(tmp_path / "report.txt")

        old_argv = sys.argv
        sys.argv = [
            "validate_samplesheet.py",
            "--input", str(sheet),
            "--output", out_tsv,
            "--report", report,
        ]
        try:
            vs.main()
        except SystemExit as e:
            pytest.fail(f"Unexpected exit: {e}")
        finally:
            sys.argv = old_argv

        validated = Path(out_tsv).read_text()
        assert "S001" in validated
        assert str(fwd) in validated

    def test_full_validation_missing_column(self, tmp_path):
        vs = self._import_vs()
        sheet = tmp_path / "bad.csv"
        sheet.write_text("sample_id\nS001\n")
        out_tsv = str(tmp_path / "validated.tsv")
        report  = str(tmp_path / "report.txt")

        old_argv = sys.argv
        sys.argv = [
            "validate_samplesheet.py",
            "--input", str(sheet),
            "--output", out_tsv,
            "--report", report,
        ]
        with pytest.raises(SystemExit):
            vs.main()
        sys.argv = old_argv

    def test_duplicate_sample_id(self, tmp_path):
        vs = self._import_vs()
        fwd = tmp_path / "S001_F.ab1"
        fwd.write_bytes(b"")

        csv_content = (
            f"sample_id,ab1_forward,ab1_reverse\n"
            f"S001,{fwd},\n"
            f"S001,{fwd},\n"
        )
        sheet = tmp_path / "dup.csv"
        sheet.write_text(csv_content)
        out_tsv = str(tmp_path / "validated.tsv")
        report  = str(tmp_path / "report.txt")

        old_argv = sys.argv
        sys.argv = [
            "validate_samplesheet.py",
            "--input", str(sheet),
            "--output", out_tsv,
            "--report", report,
        ]
        with pytest.raises(SystemExit):
            vs.main()
        sys.argv = old_argv

        report_text = Path(report).read_text()
        assert "Duplicate" in report_text


# ============================================================================
# Tests: genotype_assign module
# ============================================================================

class TestGenotypeAssign:
    """Tests for genotype assignment logic."""

    def _import_ga(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "genotype_assign",
            Path(__file__).parent.parent / "workflow" / "bin" / "genotype_assign.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_extract_genotype_underscore(self):
        ga = self._import_ga()
        gt, subgt = ga._extract_genotype("HBV_GENOTYPE_A_NC003977")
        assert gt == "A"

    def test_extract_genotype_jn_style(self):
        ga = self._import_ga()
        gt, subgt = ga._extract_genotype("JN642165_A")
        assert gt == "A"

    def test_extract_genotype_unknown(self):
        ga = self._import_ga()
        gt, subgt = ga._extract_genotype("random_sequence_no_gt")
        assert gt == "UNKNOWN"

    def test_assign_genotype_no_hits(self):
        ga = self._import_ga()
        gt, subgt, pct, hit, conf = ga.assign_genotype([], min_pct_id=90.0)
        assert gt == "UNKNOWN"
        assert conf == "LOW"

    def test_assign_genotype_high_confidence(self):
        ga = self._import_ga()
        hits = [
            {"sseqid": "HBV_GENOTYPE_B_AB073858", "pident": 95.0, "bitscore": 1000.0, "evalue": 1e-50},
            {"sseqid": "HBV_GENOTYPE_B_AB073859", "pident": 94.0, "bitscore": 950.0,  "evalue": 1e-48},
            {"sseqid": "HBV_GENOTYPE_C_AB014381", "pident": 85.0, "bitscore": 700.0,  "evalue": 1e-30},
        ]
        gt, subgt, pct, hit, conf = ga.assign_genotype(hits, min_pct_id=90.0)
        assert gt == "B"
        assert conf == "HIGH"

    def test_assign_genotype_below_threshold(self):
        ga = self._import_ga()
        hits = [
            {"sseqid": "HBV_GENOTYPE_A_NC003977", "pident": 75.0, "bitscore": 500.0, "evalue": 1e-20},
        ]
        gt, subgt, pct, hit, conf = ga.assign_genotype(hits, min_pct_id=90.0)
        assert conf == "LOW"

    def test_parse_blast_fmt6_valid(self, tmp_path):
        ga = self._import_ga()
        blast_output = (
            "query1\tHBV_GENOTYPE_C_AB014381\t96.5\t1000\t35\t2\t1\t1000\t1\t1000\t1e-200\t1800\n"
            "query1\tHBV_GENOTYPE_B_AB073858\t89.0\t1000\t110\t5\t1\t1000\t1\t1000\t1e-150\t1500\n"
        )
        f = tmp_path / "blast.txt"
        f.write_text(blast_output)
        hits = ga.parse_blast_fmt6(str(f))
        assert len(hits) == 2
        assert hits[0]["pident"] == 96.5
        assert hits[0]["bitscore"] == 1800.0

    def test_parse_blast_fmt6_empty(self, tmp_path):
        ga = self._import_ga()
        f = tmp_path / "empty.txt"
        f.write_text("")
        hits = ga.parse_blast_fmt6(str(f))
        assert hits == []


# ============================================================================
# Tests: trim_reads module
# ============================================================================

class TestTrimReads:
    """Tests for quality trimming logic."""

    def _import_tr(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "trim_reads",
            Path(__file__).parent.parent / "workflow" / "bin" / "trim_reads.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_find_trim_positions_no_low_quality(self):
        tr = self._import_tr()
        quals = [30] * 100
        left, right = tr.find_trim_positions(quals, fixed_trim=5, min_q=20)
        assert left == 5
        assert right == 95

    def test_find_trim_positions_low_quality_ends(self):
        tr = self._import_tr()
        # fixed_trim=20 removes first 20 and last 20 bases first → initial bounds [20, 80]
        # The 10 low-quality bases at positions 90-99 are already inside the fixed trim zone
        # so right boundary stays at 80 (first good-quality position from right boundary)
        quals = [5] * 30 + [30] * 60 + [5] * 10
        left, right = tr.find_trim_positions(quals, fixed_trim=20, min_q=20)
        assert left == 30  # advance past the 30 low-quality bases from left
        assert right == 80  # fixed trim already covers the low-quality tail

    def test_count_ambiguous(self):
        tr = self._import_tr()
        assert tr.count_ambiguous("ATCGRYSWKM") == 6
        assert tr.count_ambiguous("ATCG") == 0

    def test_full_trim_pass(self, tmp_path):
        tr = self._import_tr()
        # Create a 300-bp sequence with good quality
        seq = "ATCG" * 75  # 300 bp
        fasta = tmp_path / "read.fasta"
        fasta.write_text(f">SAMP001_forward_trimmed status=PASS\n{seq}\n")

        # Quality TSV with all Q=30
        qual_lines = ["sample_id\tdirection\tposition\tbase\tphred_quality\tis_ambiguous"]
        for i, b in enumerate(seq, start=1):
            qual_lines.append(f"SAMP001\tforward\t{i}\t{b}\t30\tfalse")
        qual_tsv = tmp_path / "quality.tsv"
        qual_tsv.write_text("\n".join(qual_lines) + "\n")

        out_fasta = str(tmp_path / "trimmed.fasta")
        out_stats = str(tmp_path / "stats.tsv")

        old_argv = sys.argv
        sys.argv = [
            "trim_reads.py",
            "--fasta", str(fasta),
            "--quality", str(qual_tsv),
            "--sample-id", "SAMP001",
            "--direction", "forward",
            "--min-quality", "20",
            "--min-length", "50",
            "--trim-ends", "10",
            "--out-fasta", out_fasta,
            "--out-stats", out_stats,
        ]
        try:
            tr.main()
        finally:
            sys.argv = old_argv

        stats = Path(out_stats).read_text()
        assert "PASS" in stats

    def test_full_trim_fail_too_short(self, tmp_path):
        tr = self._import_tr()
        # Only 30 bp – after trimming ends (20 each) → 0 bp → FAIL
        seq = "A" * 30
        fasta = tmp_path / "short.fasta"
        fasta.write_text(f">SAMP001_forward_trimmed status=PASS\n{seq}\n")

        qual_lines = ["sample_id\tdirection\tposition\tbase\tphred_quality\tis_ambiguous"]
        for i, b in enumerate(seq, start=1):
            qual_lines.append(f"SAMP001\tforward\t{i}\t{b}\t30\tfalse")
        qual_tsv = tmp_path / "quality.tsv"
        qual_tsv.write_text("\n".join(qual_lines) + "\n")

        out_fasta = str(tmp_path / "trimmed.fasta")
        out_stats = str(tmp_path / "stats.tsv")

        old_argv = sys.argv
        sys.argv = [
            "trim_reads.py",
            "--fasta", str(fasta),
            "--quality", str(qual_tsv),
            "--sample-id", "SAMP001",
            "--direction", "forward",
            "--min-quality", "20",
            "--min-length", "200",
            "--trim-ends", "20",
            "--out-fasta", out_fasta,
            "--out-stats", out_stats,
        ]
        try:
            tr.main()
        finally:
            sys.argv = old_argv

        stats = Path(out_stats).read_text()
        assert "FAIL" in stats


# ============================================================================
# Tests: annotate_variants module
# ============================================================================

class TestAnnotateVariants:
    """Tests for variant annotation logic."""

    def _import_av(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "annotate_variants",
            Path(__file__).parent.parent / "workflow" / "bin" / "annotate_variants.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_domain_label(self):
        av = self._import_av()
        assert av._domain_label(600) == "S"
        assert av._domain_label(300) == "preS1"
        assert av._domain_label(500) == "preS2"
        assert av._domain_label(900) == "MHR"
        assert av._domain_label(940) == "a_determinant"

    def test_translate_frame_simple(self):
        av = self._import_av()
        # ATG = Met
        result = av.translate_frame("ATGAAATTT")
        assert result[0] == "M"
        assert result[1] == "K"
        assert result[2] == "F"

    def test_annotate_no_variants(self):
        av = self._import_av()
        # identical sequences → no variants
        seq = "A" * 1250
        variants = av.annotate("S001", seq, seq, catalogue=[])
        assert variants == []

    def test_annotate_g145r_detection(self):
        av = self._import_av()
        # Build catalogue entry for G145R
        catalogue = [
            {
                "aa_position": 145,
                "ref_aa": "G",
                "alt_aa": "R",
                "region": "a_determinant",
                "mechanism": "immune_escape",
                "evidence_level": "A",
                "notes": "G145R test",
            }
        ]
        # Create ref and consensus sequences
        # The S gene starts at position 550 in our domain map
        # G is at codon 145 of S → nt position in S = (145-1)*3 = 432
        # → nt position in amplicon = 550 + 432 = 982
        ref_seq = "A" * 1250
        cons_seq = list("A" * 1250)
        # Place GGG (Gly) in reference S ORF position 145
        # Place CGG (Arg) in consensus
        s_codon_start = 550 + (145 - 1) * 3  # = 982
        # ref has AAA there → translated will be Lys unless we set it up properly
        # For simplicity, build partial sequence that differs at one codon
        # The translate function works on the full S region (nt 550 onwards)
        # We'll test annotate() with minimal sequences to detect the diff
        # Use a simple 30-aa test: just ensure annotation is triggered
        ref_s_region = "GGG" * 150  # all Gly codons
        cons_s_region = "GGG" * 144 + "CGG" + "GGG" * 5  # Arg at position 145

        ref_full  = "A" * 550 + ref_s_region  + "A" * (1250 - 550 - len(ref_s_region))
        cons_full = "A" * 550 + cons_s_region + "A" * (1250 - 550 - len(cons_s_region))

        variants = av.annotate("S001", cons_full, ref_full, catalogue)
        positions = [v["aa_position"] for v in variants]
        assert 145 in positions
        match = next(v for v in variants if v["aa_position"] == 145)
        assert match["evidence_level"] == "A"
        assert match["catalogued"] == "yes"


# ============================================================================
# Tests: generate_report module
# ============================================================================

class TestGenerateReport:
    """Tests for report generation."""

    def _import_gr(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "generate_report",
            Path(__file__).parent.parent / "workflow" / "bin" / "generate_report.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_generate_tsv(self):
        gr = self._import_gr()
        content = gr.generate_tsv(
            sample_id="S001",
            genotype_rows=[{"genotype": "B", "subgenotype": "B2", "pct_identity": "95.0",
                            "blast_hit": "ref", "confidence": "HIGH", "note": ""}],
            variant_rows=[],
            trim_rows=[],
            consensus_rows=[],
            generated_at="2026-01-01 00:00:00 UTC",
        )
        assert "S001" in content
        assert "GENOTYPE" in content
        assert "B" in content

    def test_generate_markdown(self):
        gr = self._import_gr()
        content = gr.generate_markdown(
            sample_id="S001",
            cons_seq="ATCGATCG",
            genotype_rows=[{"genotype": "C", "subgenotype": "", "confidence": "HIGH",
                            "pct_identity": "97.0", "blast_hit": "ref", "note": ""}],
            variant_rows=[],
            trim_rows=[],
            consensus_rows=[],
            generated_at="2026-01-01 00:00:00 UTC",
        )
        # [EN] Header is now bilingual VI/EN
        # [VI] Tiêu đề song ngữ VI/EN
        assert "HBV-OBI Sanger Analysis Report" in content
        assert "S001" in content
        assert "ATCGATCG" in content

    def test_generate_html(self):
        gr = self._import_gr()
        content = gr.generate_html(
            sample_id="S001",
            cons_header="S001 method=paired",
            cons_seq="ATCGATCG",
            genotype_rows=[{"genotype": "A", "subgenotype": "A2", "pct_identity": "98.0",
                            "blast_hit": "NC_003977", "confidence": "HIGH", "note": ""}],
            variant_rows=[{"notation": "sG145R", "domain": "a_determinant",
                           "mechanism": "immune_escape", "evidence_level": "A",
                           "catalogued": "yes", "notes": "G145R test"}],
            trim_rows=[],
            consensus_rows=[],
            generated_at="2026-01-01 00:00:00 UTC",
        )
        assert "<!DOCTYPE html>" in content
        assert "S001" in content
        assert "sG145R" in content
        assert "Review recommended" in content

    def test_full_report_generation_html(self, tmp_path):
        gr = self._import_gr()
        # Create input files
        cons = tmp_path / "consensus.fasta"
        cons.write_text(">S001 method=single_read\nATCGATCG\n")

        gt = tmp_path / "genotype.tsv"
        gt.write_text(
            "sample_id\tgenotype\tsubgenotype\tpct_identity\tblast_hit\tconfidence\tnote\n"
            "S001\tB\tB2\t95.0\tref_B\tHIGH\t\n"
        )
        var = tmp_path / "variants.tsv"
        var.write_text(
            "sample_id\taa_position\tref_aa\talt_aa\tnotation\tdomain\t"
            "mechanism\tevidence_level\tnotes\tcatalogued\n"
        )
        trim = tmp_path / "trim.tsv"
        trim.write_text(
            "sample_id\tdirection\traw_length\ttrimmed_length\ttrim_left\ttrim_right\t"
            "mean_quality\tambiguous_bases\tqc_status\n"
            "S001\tforward\t500\t450\t25\t475\t28.5\t3\tPASS\n"
        )
        cs = tmp_path / "consensus_stats.tsv"
        cs.write_text(
            "sample_id\tmethod\tconsensus_length\toverlap_length\t"
            "conflict_positions\tconflict_policy\tn_pass_reads\n"
            "S001\tsingle_read_forward\t450\t0\t0\tiupac\t1\n"
        )
        out = tmp_path / "report.html"

        old_argv = sys.argv
        sys.argv = [
            "generate_report.py",
            "--sample-id", "S001",
            "--consensus", str(cons),
            "--genotype", str(gt),
            "--variants", str(var),
            "--trim-stats", str(trim),
            "--consensus-stats", str(cs),
            "--format", "html",
            "--out", str(out),
        ]
        try:
            gr.main()
        finally:
            sys.argv = old_argv

        result = out.read_text()
        assert "<!DOCTYPE html>" in result
        assert "S001" in result
        assert "Genotype B" in result or "B" in result
