import difflib
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


def compare_files(file1, file2, verbose=False):
    """
    Compare generated and standard files while ignoring meta-header lines
    beginning with ## and normalizing path components.

    The #CHROM line and all variant records are compared exactly after
    path normalization.
    """

    def normalize_line(line):
        import re

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
        with open(filename, "r", encoding="utf-8-sig") as f:
            return [
                normalize_line(line)
                for line in f
                if not line.startswith("##")
            ]

    content1 = read_file(file1)
    content2 = read_file(file2)

    if len(content1) != len(content2):
        if verbose:
            print(
                "Files have different number of lines: "
                f"{len(content1)} vs {len(content2)}"
            )

            print("\nContent of file1:")
            for line in content1:
                print(repr(line))

            print("\nContent of file2:")
            for line in content2:
                print(repr(line))

        return False

    differences = []

    for i, (line1, line2) in enumerate(zip(content1, content2), 1):
        if line1 != line2:
            if verbose:
                print(f"\nDifference at line {i}:")
                print(f"File 1: {repr(line1)}")
                print(f"File 2: {repr(line2)}")

                if len(line1) != len(line2):
                    print(
                        "Line lengths differ: "
                        f"{len(line1)} vs {len(line2)}"
                    )

                for j, s in enumerate(difflib.ndiff(line1, line2)):
                    if s[0] == " ":
                        continue
                    if s[0] == "-":
                        print(
                            f"Delete '{s[-1]}' from position {j}"
                        )
                    elif s[0] == "+":
                        print(
                            f"Add '{s[-1]}' to position {j}"
                        )

            differences.append((i, line1, line2))

    return len(differences) == 0


class TestOctopusV:
    def setup_method(self):
        """Runs before each test; ensure the output directory exists."""
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    def _prepare_corrected_merge_inputs(self):
        """
        Generate the three corrected SVCF files used by merge tests.

        The fixed order is:
            sniffles, svim, pbsv

        Merge output provenance ordering depends on input order, so all merge
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

        This protects sample-column ordering and source/evidence provenance.
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