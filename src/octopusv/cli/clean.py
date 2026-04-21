"""octopusv clean: fix broken VCFs so Truvari and bcftools can read them.

This command sanitizes problematic INFO fields, fills in missing header
definitions, fixes invalid genotypes/SVLEN, optionally harmonizes chromosome
naming against a reference FASTA, then sorts, bgzips and tabix-indexes the
output. It never filters variants — all records are preserved.
"""

import gzip
import logging
import re
import shutil
import subprocess
from pathlib import Path

import typer

logger = logging.getLogger(__name__)

# Standard chromosomes for naming harmonization (main contigs only).
STD_AUTOSOMES = {str(i) for i in range(1, 23)}
STD_SEX = {"X", "Y"}
STD_MITO = {"M", "MT"}
STD_CORE = STD_AUTOSOMES | STD_SEX  # mito handled separately (M vs MT)

# Regex to strip RNAMES=...; from INFO field. RNAMES values often contain
# commas, which break strict VCF parsing in bcftools and Truvari.
_RNAMES_PATTERN = re.compile(r";?RNAMES=.*?(?=;[A-Z_][A-Z_0-9]*=|$)")

# Regex for header definition lines.
_INFO_ID_PATTERN = re.compile(r"##INFO=<ID=([^,]+)")
_FORMAT_ID_PATTERN = re.compile(r"##FORMAT=<ID=([^,]+)")
_CONTIG_ID_PATTERN = re.compile(r"##contig=<ID=([^,>]+)")


# --------------------------------------------------------------------------- #
# Dependency check
# --------------------------------------------------------------------------- #
def _check_external_tools():
    """Ensure bcftools, bgzip and tabix are available on PATH."""
    missing = [t for t in ("bcftools", "bgzip", "tabix") if shutil.which(t) is None]
    if missing:
        typer.echo(
            f"Error: required external tool(s) not found in PATH: {', '.join(missing)}",
            err=True,
        )
        typer.echo(
            "Install via: conda install -c bioconda bcftools htslib",
            err=True,
        )
        raise typer.Exit(code=1)


# --------------------------------------------------------------------------- #
# FASTA + chromosome harmonization (Step 2)
# --------------------------------------------------------------------------- #
def _read_fasta_contigs(fasta_path: Path) -> set[str]:
    """Read contig names from a FASTA file by scanning '>' header lines only."""
    contigs = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # '>chr1 any extra description' -> 'chr1'
                contigs.add(line[1:].split()[0])
    return contigs


def _detect_fasta_style(contigs: set[str]) -> str:
    """Return 'chr' if FASTA uses chr-prefixed names, 'nochr' otherwise."""
    chr_hits = sum(1 for c in ("chr1", "chr2", "chr3") if c in contigs)
    nochr_hits = sum(1 for c in ("1", "2", "3") if c in contigs)
    return "chr" if chr_hits >= nochr_hits else "nochr"


def _harmonize_chrom(chrom: str, style: str, fasta_contigs: set[str]) -> str:
    """Convert a chromosome name to match FASTA style.

    Only main chromosomes (1-22, X, Y, M/MT) are modified; decoys, alts,
    and unplaced contigs are left untouched. If the harmonized name is
    not present in the FASTA, return the original name unchanged.
    """
    stripped = chrom[3:] if chrom.startswith("chr") else chrom

    # Mitochondrial: pick whichever form (M or MT) the FASTA actually has.
    if stripped in STD_MITO:
        if style == "chr":
            for cand in ("chrM", "chrMT"):
                if cand in fasta_contigs:
                    return cand
        else:
            for cand in ("MT", "M"):
                if cand in fasta_contigs:
                    return cand
        return chrom

    # Autosomes + sex chromosomes.
    if stripped in STD_CORE:
        target = f"chr{stripped}" if style == "chr" else stripped
        if target in fasta_contigs:
            return target
        return chrom

    # Decoys / alts / unplaced: leave as-is.
    return chrom


def _process_header_contig_line(line: str, style: str, fasta_contigs: set[str]) -> str:
    """Rewrite the ID in ##contig=<ID=...> header lines."""
    m = _CONTIG_ID_PATTERN.match(line)
    if not m:
        return line
    old = m.group(1)
    new = _harmonize_chrom(old, style, fasta_contigs)
    if new != old:
        return line.replace(f"ID={old}", f"ID={new}", 1)
    return line


def _harmonize_info_chroms(info: str, style: str, fasta_contigs: set[str]) -> str:
    """Harmonize chromosome values inside INFO fields like CHR2=/CHROM2=."""
    parts = []
    for kv in info.split(";"):
        if "=" in kv:
            key, value = kv.split("=", 1)
            if key in ("CHR2", "CHROM2"):
                value = _harmonize_chrom(value, style, fasta_contigs)
            parts.append(f"{key}={value}")
        else:
            parts.append(kv)
    return ";".join(parts)


# --------------------------------------------------------------------------- #
# Header processing (Step 3)
# --------------------------------------------------------------------------- #
def _is_valid_info_id(info_id: str) -> bool:
    """Check whether an INFO/FORMAT ID is free of VCF-breaking characters."""
    return not (
        "," in info_id
        or " " in info_id
        or "<" in info_id
        or ">" in info_id
    )


def _extract_defined_fields(header_lines: list[str]) -> tuple[set[str], set[str]]:
    """Collect IDs of INFO and FORMAT definitions already present in header."""
    info_ids = set()
    format_ids = set()
    for line in header_lines:
        m = _INFO_ID_PATTERN.match(line)
        if m:
            info_ids.add(m.group(1))
            continue
        m = _FORMAT_ID_PATTERN.match(line)
        if m:
            format_ids.add(m.group(1))
    return info_ids, format_ids


def _find_undefined_fields(
    data_lines: list[str],
    defined_info: set[str],
    defined_format: set[str],
) -> tuple[set[str], set[str]]:
    """Scan data lines for INFO/FORMAT keys not declared in header.

    Invalid-character keys are skipped here (they will be sanitized during
    data cleaning, and the sanitized form will fall back to Number=.,Type=String
    auto-stubs if still unknown).
    """
    undefined_info = set()
    undefined_format = set()

    for line in data_lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        # INFO column (index 7).
        info = fields[7]
        if info and info != ".":
            for item in info.split(";"):
                key = item.split("=", 1)[0] if "=" in item else item
                if key and key not in defined_info and _is_valid_info_id(key):
                    undefined_info.add(key)

        # FORMAT column (index 8).
        if len(fields) > 8:
            for key in fields[8].split(":"):
                if key and key not in defined_format and _is_valid_info_id(key):
                    undefined_format.add(key)

    return undefined_info, undefined_format


def _fix_headers(
    header_lines: list[str],
    undefined_info: set[str],
    undefined_format: set[str],
    will_add_gt: bool,
    will_add_svlen: bool,
) -> list[str]:
    """Fix header lines.

    Actions:
        - drop RNAMES INFO definitions
        - normalize PL to Number=G,Type=Integer
        - normalize GT to Number=1,Type=String (add if missing)
        - normalize SVLEN to Number=1,Type=Integer (add if missing)
        - add Number=.,Type=String stubs for other undefined INFO/FORMAT IDs

    The last element of header_lines must be the #CHROM column header.
    """
    fixed = []
    seen_gt = False
    seen_svlen = False
    seen_pl = False

    for line in header_lines:
        # Drop RNAMES INFO definition (field is also stripped from data).
        m = _INFO_ID_PATTERN.match(line)
        if m and m.group(1) == "RNAMES":
            continue

        if line.startswith("##FORMAT=<ID=PL,"):
            line = (
                '##FORMAT=<ID=PL,Number=G,Type=Integer,'
                'Description="Normalized, Phred-scaled likelihoods for genotypes">'
            )
            seen_pl = True

        if line.startswith("##FORMAT=<ID=GT,"):
            if seen_gt:
                # Duplicate GT; skip to keep exactly one canonical definition.
                continue
            line = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
            seen_gt = True

        if line.startswith("##INFO=<ID=SVLEN,"):
            if seen_svlen:
                continue
            line = '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">'
            seen_svlen = True

        fixed.append(line)

    # Additional header lines to insert before #CHROM.
    additions = []

    if not seen_gt and (will_add_gt or "GT" in undefined_format):
        additions.append(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        )
    if not seen_svlen and (will_add_svlen or "SVLEN" in undefined_info):
        additions.append(
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">'
        )
    if not seen_pl and "PL" in undefined_format:
        additions.append(
            '##FORMAT=<ID=PL,Number=G,Type=Integer,'
            'Description="Normalized, Phred-scaled likelihoods for genotypes">'
        )

    # Auto-stubs for other undefined fields. GT/SVLEN/PL/RNAMES handled above.
    for key in sorted(undefined_info):
        if key in ("SVLEN", "RNAMES"):
            continue
        additions.append(
            f'##INFO=<ID={key},Number=.,Type=String,'
            f'Description="Auto-generated definition for {key}">'
        )
    for key in sorted(undefined_format):
        if key in ("GT", "PL"):
            continue
        additions.append(
            f'##FORMAT=<ID={key},Number=.,Type=String,'
            f'Description="Auto-generated definition for {key}">'
        )

    # Insert additions just before the final #CHROM header line.
    if fixed and fixed[-1].startswith("#CHROM"):
        fixed = fixed[:-1] + additions + [fixed[-1]]
    else:
        fixed.extend(additions)

    return fixed


# --------------------------------------------------------------------------- #
# Data line processing (Step 4)
# --------------------------------------------------------------------------- #
def _sanitize_info_value(value: str) -> str:
    """Remove characters that break INFO parsing (semicolons, whitespace)."""
    value = value.replace(";", "_")
    value = re.sub(r"\s+", "_", value)
    return value


def _clean_info_field(info: str) -> str:
    """Remove RNAMES= entirely and sanitize other problematic keys/values."""
    info = _RNAMES_PATTERN.sub("", info)
    if info.startswith(";"):
        info = info[1:]
    if not info:
        return "."

    parts = []
    for kv in info.split(";"):
        if "=" in kv:
            key, value = kv.split("=", 1)
            if not _is_valid_info_id(key):
                key = re.sub(r"[,\s<>]", "_", key)
                value = _sanitize_info_value(value)
            elif ";" in value or re.search(r"\s", value):
                value = _sanitize_info_value(value)
            parts.append(f"{key}={value}")
        else:
            # Flag-style entry (e.g. IMPRECISE); keep as-is.
            parts.append(kv)
    return ";".join(parts) if parts else "."


def _data_needs_gt(data_lines: list[str]) -> bool:
    """Check whether any record will have GT inserted into FORMAT."""
    for line in data_lines:
        fields = line.split("\t")
        if len(fields) >= 9 and "GT" not in fields[8].split(":"):
            return True
    return False


def _data_needs_svlen(data_lines: list[str]) -> bool:
    """Check whether any record will have SVLEN filled in."""
    for line in data_lines:
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        info = fields[7]
        svlen_found = False
        for kv in info.split(";"):
            if kv.startswith("SVLEN="):
                v = kv.split("=", 1)[1]
                if v and v != ".":
                    svlen_found = True
                break
        if not svlen_found:
            return True
    return False


def _fix_svlen(info: str, pos: str) -> str:
    """Fill in SVLEN when missing or '.', using SVTYPE and END if available."""
    info_dict = {}
    flags = []
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            info_dict[k] = v
        elif kv:
            flags.append(kv)

    svlen = info_dict.get("SVLEN")
    svtype = info_dict.get("SVTYPE")
    end_str = info_dict.get("END")

    if svlen is None or svlen == ".":
        try:
            curr_pos = int(pos)
        except ValueError:
            return info

        if end_str is not None and end_str.isdigit():
            end = int(end_str)
            if svtype == "DEL":
                new_svlen = -(end - curr_pos + 1)
            elif svtype == "INS":
                new_svlen = max(1, end - curr_pos)
            elif svtype in ("DUP", "INV"):
                new_svlen = max(1, end - curr_pos + 1)
            else:
                new_svlen = 0
        else:
            if svtype == "INS":
                new_svlen = 1
            elif svtype == "DEL":
                new_svlen = -1
            else:
                new_svlen = 0

        info_dict["SVLEN"] = str(new_svlen)

    pieces = [f"{k}={v}" for k, v in info_dict.items()] + flags
    return ";".join(pieces) if pieces else "."


def _fix_genotype(format_str: str, sample_str: str) -> tuple[str, str]:
    """Ensure GT is first in FORMAT and has a valid value."""
    fmt = format_str.split(":")
    smp = sample_str.split(":")

    if "GT" not in fmt:
        fmt.insert(0, "GT")
        smp.insert(0, "1/1")
    elif fmt[0] != "GT":
        idx = fmt.index("GT")
        fmt.insert(0, fmt.pop(idx))
        if idx < len(smp):
            smp.insert(0, smp.pop(idx))
        else:
            smp.insert(0, "1/1")

    # Replace invalid GT calls with a benign homozygous-alt call.
    if smp[0] in ("./.", ".", "./1", "1/.", ".|.", ".|1", "1|."):
        smp[0] = "1/1"

    # Pad sample if shorter than format.
    while len(smp) < len(fmt):
        smp.append(".")

    return ":".join(fmt), ":".join(smp)


# --------------------------------------------------------------------------- #
# Main pipeline
# --------------------------------------------------------------------------- #
def _process_vcf(
    input_vcf: Path,
    output_path: Path,
    fasta_path: Path | None,
):
    """Run the full clean pipeline, producing a bgzipped + tabix-indexed VCF."""
    # Step 2 setup: load FASTA contigs if a reference was provided.
    fasta_contigs: set[str] = set()
    style: str | None = None
    if fasta_path is not None:
        typer.echo(f"Reading contigs from reference FASTA: {fasta_path}")
        fasta_contigs = _read_fasta_contigs(fasta_path)
        if not fasta_contigs:
            typer.echo(f"Error: no contigs found in FASTA {fasta_path}", err=True)
            raise typer.Exit(code=1)
        style = _detect_fasta_style(fasta_contigs)
        typer.echo(f"Detected FASTA naming style: {style}")

    # Read input VCF (plain or gzip).
    opener = gzip.open if str(input_vcf).endswith(".gz") else open
    header_lines: list[str] = []
    data_lines: list[str] = []
    with opener(input_vcf, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                # Step 2a: harmonize ##contig=<ID=...> header lines.
                if fasta_path is not None and line.startswith("##contig=<ID="):
                    line = _process_header_contig_line(line, style, fasta_contigs)
                header_lines.append(line)
            else:
                data_lines.append(line)

    if not header_lines or not header_lines[-1].startswith("#CHROM"):
        typer.echo("Error: input VCF is missing a #CHROM header line.", err=True)
        raise typer.Exit(code=1)

    typer.echo(f"Read {len(data_lines)} variant records")

    # Step 2b: harmonize CHROM column + INFO CHR2/CHROM2 on data lines.
    if fasta_path is not None:
        harmonized = []
        for line in data_lines:
            fields = line.split("\t")
            if len(fields) < 8:
                harmonized.append(line)
                continue
            fields[0] = _harmonize_chrom(fields[0], style, fasta_contigs)
            fields[7] = _harmonize_info_chroms(fields[7], style, fasta_contigs)
            harmonized.append("\t".join(fields))
        data_lines = harmonized

    # Step 3: fix headers (drop RNAMES, auto-stub undefined fields, normalize
    # GT/SVLEN/PL). We precompute whether data will induce GT/SVLEN additions
    # so the header stays consistent with the data we are about to write.
    defined_info, defined_format = _extract_defined_fields(header_lines)
    undefined_info, undefined_format = _find_undefined_fields(
        data_lines, defined_info, defined_format
    )
    will_add_gt = _data_needs_gt(data_lines)
    will_add_svlen = _data_needs_svlen(data_lines)
    fixed_headers = _fix_headers(
        header_lines, undefined_info, undefined_format, will_add_gt, will_add_svlen
    )

    # Step 4: fix each data line.
    cleaned_data: list[str] = []
    for line in data_lines:
        fields = line.split("\t")
        if len(fields) < 8:
            cleaned_data.append(line)  # malformed row; keep verbatim
            continue

        # FILTER -> PASS so Truvari does not drop it by default.
        fields[6] = "PASS"

        # INFO: strip RNAMES and sanitize illegal characters.
        fields[7] = _clean_info_field(fields[7])

        # INFO: fill missing SVLEN based on SVTYPE + END.
        fields[7] = _fix_svlen(fields[7], fields[1])

        # FORMAT + first SAMPLE column: ensure a valid GT.
        if len(fields) >= 10:
            fields[8], fields[9] = _fix_genotype(fields[8], fields[9])

        cleaned_data.append("\t".join(fields))

    # Determine output prefix (tolerate .vcf, .vcf.gz, or no extension).
    out_str = str(output_path)
    if out_str.endswith(".vcf.gz"):
        output_prefix = out_str[:-7]
    elif out_str.endswith(".vcf"):
        output_prefix = out_str[:-4]
    else:
        output_prefix = out_str

    temp_vcf = Path(f"{output_prefix}_temp_unsorted.vcf")
    sorted_vcf = Path(f"{output_prefix}_sorted.vcf")
    final_gz = Path(f"{output_prefix}.vcf.gz")

    # Write intermediate plain VCF.
    with open(temp_vcf, "w") as f:
        for line in fixed_headers:
            f.write(line + "\n")
        for line in cleaned_data:
            f.write(line + "\n")

    # Sort using bcftools (respects contig order from header).
    typer.echo("Sorting with bcftools...")
    try:
        subprocess.run(
            ["bcftools", "sort", str(temp_vcf), "-o", str(sorted_vcf)],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        typer.echo("Error during bcftools sort:", err=True)
        typer.echo(e.stderr, err=True)
        raise typer.Exit(code=1) from e

    # bgzip compress + tabix index.
    typer.echo("Compressing and indexing...")
    try:
        with open(final_gz, "wb") as out:
            subprocess.run(["bgzip", "-c", str(sorted_vcf)], check=True, stdout=out)
        subprocess.run(
            ["tabix", "-p", "vcf", str(final_gz)],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        typer.echo(f"Error during compression/indexing: {e}", err=True)
        raise typer.Exit(code=1) from e

    # Cleanup intermediates.
    temp_vcf.unlink(missing_ok=True)
    sorted_vcf.unlink(missing_ok=True)

    typer.echo(f"Done. Output: {final_gz}")
    typer.echo(f"Index:  {final_gz}.tbi")


# --------------------------------------------------------------------------- #
# Typer command
# --------------------------------------------------------------------------- #
def clean(
    input_vcf: Path | None = typer.Argument(
        None,
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="Input VCF file to clean (positional).",
    ),
    output: Path | None = typer.Argument(
        None,
        dir_okay=False,
        resolve_path=True,
        help="Output file path (positional). .vcf.gz will be produced.",
    ),
    input_option: Path | None = typer.Option(
        None,
        "--input-file",
        "-i",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="Input VCF file to clean.",
    ),
    output_option: Path | None = typer.Option(
        None,
        "--output-file",
        "-o",
        dir_okay=False,
        resolve_path=True,
        help="Output file path. .vcf.gz will be produced.",
    ),
    genome: Path | None = typer.Option(
        None,
        "--genome",
        "-g",
        exists=True,
        dir_okay=False,
        resolve_path=True,
        help="Reference genome FASTA (absolute path). If provided, chromosome "
             "naming in the VCF will be harmonized against this FASTA.",
    ),
):
    """Clean a broken VCF so that bcftools/Truvari can parse it.

    Removes RNAMES, sanitizes illegal INFO characters, fills missing GT and
    SVLEN, adds missing header definitions, optionally harmonizes chromosome
    names against a reference FASTA, then sorts, bgzips and tabix-indexes
    the output. No variants are filtered.
    """
    # Resolve input.
    if input_vcf and input_option:
        typer.echo(
            "Error: specify input either as positional argument or -i, not both.",
            err=True,
        )
        raise typer.Exit(code=1)
    input_file = input_vcf or input_option
    if input_file is None:
        typer.echo("Error: input file is required.", err=True)
        raise typer.Exit(code=1)

    # Resolve output.
    if output and output_option:
        typer.echo(
            "Error: specify output either as positional argument or -o, not both.",
            err=True,
        )
        raise typer.Exit(code=1)
    output_file = output or output_option
    if output_file is None:
        typer.echo("Error: output file is required.", err=True)
        raise typer.Exit(code=1)

    _check_external_tools()
    _process_vcf(input_file, output_file, genome)