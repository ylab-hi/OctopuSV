import subprocess
import pytest
import os
import difflib
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Get absolute path of project root directory
ROOT_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

# Define directory structure using absolute paths
TEST_DATA_DIR = os.path.join(ROOT_DIR, "tests", "data")
INPUT_DIR = os.path.join(TEST_DATA_DIR, "input")
STANDARD_DIR = os.path.join(TEST_DATA_DIR, "standard")
OUTPUT_DIR = os.path.join(TEST_DATA_DIR, "output")


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
            cwd=ROOT_DIR
        )
        if verbose and result.stdout:
            logger.info(f"Command stdout:\n{result.stdout}")
        if result.stderr:
            logger.warning(f"Command stderr:\n{result.stderr}")
        result.check_returncode()
    except subprocess.CalledProcessError as e:
        logger.error(f"Command execution failed with exit code {e.returncode}")
        logger.error(f"Error output:\n{e.stderr}")
        raise


def compare_files(file1, file2, verbose=False):
    """
    Enhanced file comparison function that handles various edge cases
    and provides detailed difference information.
    """

    def normalize_line(line):
        # Remove common path patterns
        import re
        # Remove absolute paths
        line = re.sub(r'/[^,\s;]*/([^/,\s;]+)', r'\1', line)
        # Remove relative paths
        line = re.sub(r'\.\.?/[^,\s;]+/([^/,\s;]+)', r'\1', line)
        return line.rstrip('\r\n')

    def read_file(filename):
        with open(filename, 'r', encoding='utf-8-sig') as f:  # handles BOM
            return [normalize_line(line) for line in f
                    if not line.startswith('##')]

    content1 = read_file(file1)
    content2 = read_file(file2)

    if len(content1) != len(content2):
        if verbose:
            print(f"Files have different number of lines: {len(content1)} vs {len(content2)}")
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

                # Show detailed character comparison
                if len(line1) != len(line2):
                    print(f"Line lengths differ: {len(line1)} vs {len(line2)}")

                # Use difflib to show exact differences
                for j, s in enumerate(difflib.ndiff(line1, line2)):
                    if s[0] == ' ':
                        continue
                    elif s[0] == '-':
                        print(f"Delete '{s[-1]}' from position {j}")
                    elif s[0] == '+':
                        print(f"Add '{s[-1]}' to position {j}")

            differences.append((i, line1, line2))

    return len(differences) == 0


class TestOctopusV:
    def setup_method(self):
        """Runs before each test; ensure the output directory exists."""
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---- three independent correct tests -----------------------------------

    def test_correct_sniffles(self):
        self._check_correct("sniffles")

    def test_correct_svim(self):
        self._check_correct("svim")

    def test_correct_pbsv(self):
        self._check_correct("pbsv")

    def _check_correct(self, caller):
        """Run `octopusv correct` on one caller VCF and compare to standard."""
        input_vcf = os.path.join(INPUT_DIR, f"{caller}.vcf")
        output_svcf = os.path.join(OUTPUT_DIR, f"{caller}.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, f"{caller}.svcf")

        assert os.path.exists(input_vcf), f"Input not found: {input_vcf}"
        assert os.path.exists(standard_svcf), f"Standard not found: {standard_svcf}"

        run_octopusv("correct", "-i", input_vcf, "-o", output_svcf, verbose=False)
        assert os.path.exists(output_svcf), f"Output not created: {output_svcf}"
        assert compare_files(output_svcf, standard_svcf, verbose=True), \
            f"{caller} correct output does not match standard"

    # ---- merge over the three corrected SVCFs ------------------------------

    def test_merge_min_support(self):
        """Merge the three corrected SVCFs with --min-support 2 (fixed order).

        The correct step is run here so the test is self-contained and a
        regression in correct also fails this test. Inputs deliberately include
        underscore-containing contigs (e.g. chrY_KI270740v1_random, NC_007605)
        to guard the CO coordinate parser against the unpack regression.
        """
        for caller in ("sniffles", "svim", "pbsv"):
            run_octopusv("correct",
                         "-i", os.path.join(INPUT_DIR, f"{caller}.vcf"),
                         "-o", os.path.join(OUTPUT_DIR, f"{caller}.svcf"),
                         verbose=False)

        output_svcf = os.path.join(OUTPUT_DIR, "min2.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "min2.svcf")
        assert os.path.exists(standard_svcf), f"Standard not found: {standard_svcf}"

        # Fixed input order: sniffles, svim, pbsv (order affects SOURCES).
        run_octopusv("merge",
                     "-i",
                     os.path.join(OUTPUT_DIR, "sniffles.svcf"),
                     os.path.join(OUTPUT_DIR, "svim.svcf"),
                     os.path.join(OUTPUT_DIR, "pbsv.svcf"),
                     "--min-support", "2",
                     "-o", output_svcf, verbose=True)

        assert os.path.exists(output_svcf), f"Output not created: {output_svcf}"
        assert compare_files(output_svcf, standard_svcf, verbose=True), \
            "Merge --min-support output does not match standard"

    # ---- svcf2vcf on the merged result -------------------------------------

    def test_svcf2vcf(self):
        """Convert the standard merged SVCF to VCF and compare."""
        input_svcf = os.path.join(STANDARD_DIR, "min2.svcf")
        output_vcf = os.path.join(OUTPUT_DIR, "min2.vcf")
        standard_vcf = os.path.join(STANDARD_DIR, "min2.vcf")

        assert os.path.exists(input_svcf), f"Input not found: {input_svcf}"
        assert os.path.exists(standard_vcf), f"Standard not found: {standard_vcf}"

        run_octopusv("svcf2vcf", "-i", input_svcf, "-o", output_vcf, verbose=False)
        assert os.path.exists(output_vcf), f"Output not created: {output_vcf}"
        assert compare_files(output_vcf, standard_vcf, verbose=True), \
            "svcf2vcf output does not match standard"

# test