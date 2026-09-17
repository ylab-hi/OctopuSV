import logging
import os
import subprocess


# Configure logging
logging.basicConfig(level=logging.DEBUG, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# Get absolute path of project root directory
ROOT_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

# Define directory structure using absolute paths
TEST_DATA_DIR = os.path.join(ROOT_DIR, "tests", "data")
INPUT_DIR = os.path.join(TEST_DATA_DIR, "input")
STANDARD_DIR = os.path.join(TEST_DATA_DIR, "standard")
OUTPUT_DIR = os.path.join(TEST_DATA_DIR, "output")


# Fixed merge input order.
# Order matters for SOURCES, SOURCE_IDS, and evidence/sample columns.
MERGE_CALLERS = ("sniffles", "svim", "pbsv")


def run_octopusv(command, *args, verbose=False):
    """
    Run octopusv command and capture output.
    If command fails, print detailed error information and raise exception.
    """
    cmd = ["octopusv", command] + list(args)

    if verbose:
        logger.info(f"Executing command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=ROOT_DIR,
        )

        if verbose and result.stdout:
            logger.info(f"Command stdout:\n{result.stdout}")

        if result.stderr:
            logger.warning(f"Command stderr:\n{result.stderr}")

        result.check_returncode()

    except subprocess.CalledProcessError as e:
        logger.error(
            f"Command execution failed with exit code {e.returncode}"
        )
        logger.error(f"Error output:\n{e.stderr}")
        raise


def compare_files(
    file1,
    file2,
    verbose=False,
    max_reported_diffs=10,
):
    """
    Compare generated and standard files while ignoring meta-header lines
    beginning with ## and normalizing path components.

    The #CHROM line and all variant records are compared exactly after
    path normalization.

    When differences are found, report a compact line-level summary instead
    of running a character-by-character diff on very long SVCF records.
    """
    import re

    def normalize_line(line):
        # Remove absolute paths while retaining the final basename.
        line = re.sub(
            r"/[^,\s;]*/([^/,\s;]+)",
            r"\1",
            line,
        )

        # Remove relative path prefixes while retaining the final basename.
        line = re.sub(
            r"\.\.?/[^,\s;]+/([^/,\s;]+)",
            r"\1",
            line,
        )

        return line.rstrip("\r\n")

    def read_file(filename):
        with open(
            filename,
            "r",
            encoding="utf-8-sig",
        ) as handle:
            return [
                normalize_line(line)
                for line in handle
                if not line.startswith("##")
            ]

    def parse_info(info_text):
        info = {}

        for item in info_text.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                info[key] = value
            elif item:
                info[item] = True

        return info

    def summarize_line(line):
        fields = line.split("\t")

        if not fields:
            return "<empty line>"

        if fields[0] == "#CHROM":
            preview = fields[:12]
            suffix = " ..." if len(fields) > 12 else ""
            return (
                f"#CHROM header with {len(fields)} columns: "
                f"{preview}{suffix}"
            )

        if len(fields) < 8:
            preview = line[:300]
            if len(line) > 300:
                preview += "..."
            return preview

        info = parse_info(fields[7])
        record_key = (
            fields[0],
            fields[1],
            fields[2],
        )

        summary = [
            f"record={record_key}",
            f"SOURCES={info.get('SOURCES', '.')}",
            f"SOURCE_IDS={info.get('SOURCE_IDS', '.')}",
        ]

        if len(fields) > 8:
            format_keys = fields[8].split(":")

            if "ID" in format_keys:
                id_index = format_keys.index("ID")
                evidence_ids = []

                for sample in fields[9:]:
                    values = sample.split(":")

                    if id_index < len(values):
                        evidence_ids.append(values[id_index])
                    else:
                        evidence_ids.append("<missing>")

                summary.append(f"evidence_IDs={evidence_ids}")

            summary.append(
                f"evidence_blocks={max(0, len(fields) - 9)}"
            )

        return "; ".join(summary)

    content1 = read_file(file1)
    content2 = read_file(file2)

    if len(content1) != len(content2):
        if verbose:
            print(
                "\nFiles have different numbers of non-meta lines:"
            )
            print(f"  generated: {len(content1)}")
            print(f"  standard:  {len(content2)}")

        return False

    difference_count = 0
    reported_count = 0

    for line_number, (line1, line2) in enumerate(
        zip(content1, content2),
        start=1,
    ):
        if line1 == line2:
            continue

        difference_count += 1

        if verbose and reported_count < max_reported_diffs:
            reported_count += 1

            print()
            print("=" * 80)
            print(
                f"Difference at non-meta line {line_number}"
            )

            print("Generated:")
            print("  " + summarize_line(line1))

            print("Standard:")
            print("  " + summarize_line(line2))

    if verbose and difference_count:
        print()
        print(
            f"Total differing lines: {difference_count}"
        )

        if difference_count > reported_count:
            print(
                f"Only the first {reported_count} differences were shown."
            )

    return difference_count == 0


class TestOctopusV:
    def setup_method(self):
        """Runs before each test; ensure the output directory exists."""
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    def _prepare_corrected_merge_inputs(self):
        """
        Generate the three corrected SVCF files used by merge tests.

        The fixed order is:
            sniffles, svim, pbsv

        Merge output source ordering depends on input order, so all merge
        regression tests use this same order.
        """
        corrected_files = []

        for caller in MERGE_CALLERS:
            input_vcf = os.path.join(
                INPUT_DIR,
                f"{caller}.vcf",
            )
            output_svcf = os.path.join(
                OUTPUT_DIR,
                f"{caller}.svcf",
            )

            assert os.path.exists(
                input_vcf
            ), f"Input not found: {input_vcf}"

            run_octopusv(
                "correct",
                "-i",
                input_vcf,
                "-o",
                output_svcf,
                verbose=False,
            )

            assert os.path.exists(
                output_svcf
            ), f"Output not created: {output_svcf}"

            corrected_files.append(output_svcf)

        return corrected_files

    # ---- three independent correct tests -----------------------------------

    def test_correct_sniffles(self):
        self._check_correct("sniffles")

    def test_correct_svim(self):
        self._check_correct("svim")

    def test_correct_pbsv(self):
        self._check_correct("pbsv")

    def _check_correct(self, caller):
        """Run `octopusv correct` on one caller VCF and compare to standard."""
        input_vcf = os.path.join(
            INPUT_DIR,
            f"{caller}.vcf",
        )
        output_svcf = os.path.join(
            OUTPUT_DIR,
            f"{caller}.svcf",
        )
        standard_svcf = os.path.join(
            STANDARD_DIR,
            f"{caller}.svcf",
        )

        assert os.path.exists(
            input_vcf
        ), f"Input not found: {input_vcf}"

        assert os.path.exists(
            standard_svcf
        ), f"Standard not found: {standard_svcf}"

        run_octopusv(
            "correct",
            "-i",
            input_vcf,
            "-o",
            output_svcf,
            verbose=False,
        )

        assert os.path.exists(
            output_svcf
        ), f"Output not created: {output_svcf}"

        assert compare_files(
            output_svcf,
            standard_svcf,
            verbose=True,
        ), f"{caller} correct output does not match standard"

    # ---- merge over the three corrected SVCFs ------------------------------

    def test_merge_min_support(self):
        """
        Merge the three corrected SVCFs with --min-support 2.

        Inputs deliberately include underscore-containing contigs
        (e.g. chrY_KI270740v1_random, NC_007605) to guard the CO coordinate
        parser against the unpack regression.
        """
        corrected_files = self._prepare_corrected_merge_inputs()

        output_svcf = os.path.join(
            OUTPUT_DIR,
            "min2.svcf",
        )
        standard_svcf = os.path.join(
            STANDARD_DIR,
            "min2.svcf",
        )

        assert os.path.exists(
            standard_svcf
        ), f"Standard not found: {standard_svcf}"

        run_octopusv(
            "merge",
            "-i",
            *corrected_files,
            "--min-support",
            "2",
            "-o",
            output_svcf,
            verbose=True,
        )

        assert os.path.exists(
            output_svcf
        ), f"Output not created: {output_svcf}"

        assert compare_files(
            output_svcf,
            standard_svcf,
            verbose=True,
        ), "Merge --min-support output does not match standard"

    def test_merge_union(self):
        """
        Merge all three corrected SVCFs in caller mode with --union.

        This keeps the full merged event set and therefore provides a broad
        regression check for merge grouping, representative selection,
        SOURCES, SOURCE_IDS, and caller evidence ordering.
        """
        corrected_files = self._prepare_corrected_merge_inputs()

        output_svcf = os.path.join(
            OUTPUT_DIR,
            "union.svcf",
        )
        standard_svcf = os.path.join(
            STANDARD_DIR,
            "union.svcf",
        )

        assert os.path.exists(
            standard_svcf
        ), f"Standard not found: {standard_svcf}"

        run_octopusv(
            "merge",
            "-i",
            *corrected_files,
            "--union",
            "-o",
            output_svcf,
            verbose=True,
        )

        assert os.path.exists(
            output_svcf
        ), f"Output not created: {output_svcf}"

        assert compare_files(
            output_svcf,
            standard_svcf,
            verbose=True,
        ), "Merge --union output does not match standard"

    def test_merge_sample_union(self):
        """
        Merge the same three corrected SVCFs as three independent sample-mode
        inputs with --union.

        The files are used as fixed regression inputs rather than as a
        biological sample model. Their basenames become the sample labels:
            sniffles, svim, pbsv

        This protects sample-column ordering and source/evidence mapping.
        """
        corrected_files = self._prepare_corrected_merge_inputs()

        output_svcf = os.path.join(
            OUTPUT_DIR,
            "sample_union.svcf",
        )
        standard_svcf = os.path.join(
            STANDARD_DIR,
            "sample_union.svcf",
        )

        assert os.path.exists(
            standard_svcf
        ), f"Standard not found: {standard_svcf}"

        run_octopusv(
            "merge",
            "-i",
            *corrected_files,
            "--mode",
            "sample",
            "--union",
            "-o",
            output_svcf,
            verbose=True,
        )

        assert os.path.exists(
            output_svcf
        ), f"Output not created: {output_svcf}"

        assert compare_files(
            output_svcf,
            standard_svcf,
            verbose=True,
        ), "Sample-mode merge --union output does not match standard"

    # ---- svcf2vcf on the merged result -------------------------------------

    def test_svcf2vcf(self):
        """Convert the standard merged SVCF to VCF and compare."""
        input_svcf = os.path.join(
            STANDARD_DIR,
            "min2.svcf",
        )
        output_vcf = os.path.join(
            OUTPUT_DIR,
            "min2.vcf",
        )
        standard_vcf = os.path.join(
            STANDARD_DIR,
            "min2.vcf",
        )

        assert os.path.exists(
            input_svcf
        ), f"Input not found: {input_svcf}"

        assert os.path.exists(
            standard_vcf
        ), f"Standard not found: {standard_vcf}"

        run_octopusv(
            "svcf2vcf",
            "-i",
            input_svcf,
            "-o",
            output_vcf,
            verbose=False,
        )

        assert os.path.exists(
            output_vcf
        ), f"Output not created: {output_vcf}"

        assert compare_files(
            output_vcf,
            standard_vcf,
            verbose=True,
        ), "svcf2vcf output does not match standard"
