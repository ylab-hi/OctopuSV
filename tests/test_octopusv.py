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
                for i, s in enumerate(difflib.ndiff(line1, line2)):
                    if s[0] == ' ':
                        continue
                    elif s[0] == '-':
                        print(f"Delete '{s[-1]}' from position {i}")
                    elif s[0] == '+':
                        print(f"Add '{s[-1]}' to position {i}")

            differences.append((i, line1, line2))

    return len(differences) == 0


class TestOctopusV:
    def setup_method(self):
        """Setup method runs before each test"""
        # Ensure output directory exists
        os.makedirs(OUTPUT_DIR, exist_ok=True)

    def test_correct_svim(self):
        """Test correct command with SVIM input"""
        input_vcf = os.path.join(INPUT_DIR, "svim_example.vcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_svim_example.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_svim.svcf")

        # Check if input and standard files exist
        assert os.path.exists(input_vcf), f"Input file not found: {input_vcf}"
        assert os.path.exists(standard_svcf), f"Standard file not found: {standard_svcf}"

        # Run octopusv correct command with -i and -o options
        run_octopusv("correct", "-i", input_vcf, "-o", output_svcf, verbose=False)

        # Check if output file was created
        assert os.path.exists(output_svcf), f"Output file not created: {output_svcf}"

        # Compare output file with standard file
        assert compare_files(output_svcf, standard_svcf, verbose=True), \
            "SVIM correction output does not match standard"

    def test_correct_pbsv(self):
        """Test correct command with PBSV input"""
        input_vcf = os.path.join(INPUT_DIR, "pbsv_example.vcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_pbsv_example.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_pbsv.svcf")

        # Check if input and standard files exist
        assert os.path.exists(input_vcf), f"Input file not found: {input_vcf}"
        assert os.path.exists(standard_svcf), f"Standard file not found: {standard_svcf}"

        # Run octopusv correct command with -i and -o options
        run_octopusv("correct", "-i", input_vcf, "-o", output_svcf, verbose=False)

        # Check if output file was created
        assert os.path.exists(output_svcf), f"Output file not created: {output_svcf}"

        # Compare output file with standard file
        assert compare_files(output_svcf, standard_svcf, verbose=True), \
            "PBSV correction output does not match standard"

    def test_merge_intersect(self):
        """Test merge command with intersect option"""
        input_svcf1 = os.path.join(INPUT_DIR, "pbsv.svcf")
        input_svcf2 = os.path.join(INPUT_DIR, "svim.svcf")
        output_svcf = os.path.join(OUTPUT_DIR, "out_inter.svcf")
        standard_svcf = os.path.join(STANDARD_DIR, "standard_intersection.svcf")

        # Check if input and standard files exist
        assert os.path.exists(input_svcf1), f"Input file 1 not found: {input_svcf1}"
        assert os.path.exists(input_svcf2), f"Input file 2 not found: {input_svcf2}"
        assert os.path.exists(standard_svcf), f"Standard file not found: {standard_svcf}"

        # Run octopusv merge command
        run_octopusv(
            "merge",
            "-i", input_svcf1, input_svcf2,
            "--intersect",
            "-o", output_svcf,
            verbose=True
        )

        # Check if output file was created
        assert os.path.exists(output_svcf), f"Output file not created: {output_svcf}"

        # Compare files with verbose mode
        is_same = compare_files(output_svcf, standard_svcf, verbose=True)

        # Print additional debug information if files differ
        if not is_same:
            print("\nRaw hex dump of first few lines:")
            for filename in [output_svcf, standard_svcf]:
                print(f"\n{filename}:")
                with open(filename, 'rb') as f:
                    print(f.read().hex())

        assert is_same, "Merge intersect output does not match standard"
