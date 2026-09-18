import datetime
from typing import List, Optional

from .name_mapper import NameMapper
from octopusv.utils.svcf_sample_parser import parse_svcf_sample_block
from octopusv.utils.svcf_sort import sort_events_for_output
from octopusv.utils.vcf_info import format_vcf_info_item
from octopusv.utils.svcf_schema import (
    SAMPLE_FORMAT as SAMPLE_FORMAT_V11,
    MODE_MULTI,
    mode_header,
    validate_source_id,
    validate_source_label,
    version_header,
)


class MultiSampleWriter:
    """Write SVCF output for sample mode with fixed sample-column ordering.

    In sample mode, the #CHROM trailing columns are true sample/input columns.
    Their order is fixed by the input file order or by --sample-names.

    This writer assumes each event has event.ordered_samples already prepared
    in that exact order.
    """

    def __init__(self, name_mapper: NameMapper, input_files=None):
        """Initialize writer with a NameMapper and source files for headers."""
        self.name_mapper = name_mapper
        # Header passthrough is enabled only when the orchestration layer
        # explicitly provides real input files. Direct writer use in tests or
        # libraries remains self-contained and does not assume NameMapper paths
        # are readable files.
        self.input_files = [str(path) for path in (input_files or [])]

    def write_results(self, output_file, events, contigs):
        """Write merged sample-mode SVCF results.

        Args:
            output_file: Output SVCF path.
            events: Merged events with ordered_samples attached.
            contigs: Dictionary of contig ID -> length.
        """
        ordered_events = sort_events_for_output(events, contigs.keys())

        with open(output_file, "w") as f:
            self._write_header(f, contigs)

            for event in ordered_events:
                self._write_event(f, event)

    def _write_header(self, file_handle, contigs):
        """Write SVCF header for sample mode."""
        # Extract passthrough definitions before writing any bytes so a header
        # read failure cannot leave a plausible partial header when this writer
        # is used directly outside the atomic merge CLI.
        passthrough = self._collect_passthrough_definitions()

        file_handle.write("##fileformat=VCFv4.2\n")
        file_handle.write(version_header() + "\n")
        file_handle.write(mode_header(MODE_MULTI) + "\n")

        file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
        file_handle.write(f"##fileDate={file_date}\n")
        file_handle.write("##source=OctopuSV\n")
        file_handle.write(
            "##OctopuSV_WARNING=This is SVCF format. "
            "Use 'octopusv svcf2vcf' to change back to standard VCF format "
            "before bcftools/vcftools\n"
        )

        for meta_line in passthrough.get("safe_meta_lines", []):
            file_handle.write(meta_line + "\n")

        for contig_id, contig_length in contigs.items():
            file_handle.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

        self._write_standard_definitions(file_handle)
        for category in ("alt_lines", "info_lines", "filter_lines"):
            for line in passthrough[category]:
                file_handle.write(line + "\n")

        sample_names = self.name_mapper.get_all_display_names()
        header_line = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(sample_names)
            + "\n"
        )
        file_handle.write(header_line)

    def _write_standard_definitions(self, file_handle):
        """Write core INFO and FORMAT definitions used by OctopuSV SVCF."""
        file_handle.write(
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
        )
        file_handle.write(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n'
        )
        file_handle.write(
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
        )
        file_handle.write(
            '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for end coordinate">\n'
        )
        file_handle.write(
            '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting this variant">\n'
        )
        file_handle.write(
            '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method used to detect SV">\n'
        )
        file_handle.write(
            '##INFO=<ID=RTID,Number=1,Type=String,Description="Related transcript or record ID">\n'
        )
        file_handle.write(
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n'
        )
        file_handle.write(
            '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation of the SV">\n'
        )
        file_handle.write(
            '##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names">\n'
        )
        file_handle.write(
            '##INFO=<ID=SOURCES,Number=.,Type=String,Description="List of input samples/files that support this variant">\n'
        )
        file_handle.write(
            '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="Original IDs of merged SVs from different samples/files">\n'
        )

        file_handle.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles; unavailable for synthesized sample-mode calls">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=UC,Number=1,Type=Integer,Description="Number of unique callers supporting carrier presence for this synthesized sample call">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=UV,Number=1,Type=Integer,Description="Number of unique callers contributing a valid presence vote for this synthesized sample call">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand orientation of SV">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=QV,Number=1,Type=Float,Description="Quality value">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=TY,Number=1,Type=String,Description="Type of SV">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=ID,Number=1,Type=String,Description="Unique identifier for the SV">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=SC,Number=1,Type=String,Description="Source caller/method">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=REF,Number=1,Type=String,Description="Reference allele sequence">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=ALT,Number=1,Type=String,Description="Alternate allele sequence">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinate information of the SV">\n'
        )


    def _collect_passthrough_definitions(self):
        """Preserve non-reserved INFO/FILTER/ALT definitions from inputs.

        Sample-mode records may retain representative event annotations such
        as PRECISE and caller FILTER/ALT values. Their header definitions must
        therefore remain available, while SVCF 1.1 core INFO/FORMAT semantics
        stay canonical and cannot be overridden by input headers.
        """
        from octopusv.utils.svcf_utils import (
            extract_id_from_header_line,
            extract_original_header_definitions,
            merge_safe_global_meta_lines,
        )

        canonical_info_ids = {
            "SVTYPE", "END", "SVLEN", "CHR2", "SUPPORT", "SVMETHOD",
            "RTID", "AF", "STRAND", "RNAMES", "SOURCES", "SOURCE_IDS",
        }
        collected = {
            "alt_lines": [],
            "info_lines": [],
            "filter_lines": [],
            "safe_meta_lines": [],
        }
        seen = {
            key: set()
            for key in ("alt_lines", "info_lines", "filter_lines")
        }
        safe_meta_sources = []

        category_prefix = {
            "alt_lines": "ALT",
            "info_lines": "INFO",
            "filter_lines": "FILTER",
        }

        for input_file in self.input_files:
            header = extract_original_header_definitions(input_file)
            safe_meta_sources.append(
                (str(input_file), header.get("other_lines", []))
            )
            for category, prefix in category_prefix.items():
                for line in header.get(category, []):
                    item_id = extract_id_from_header_line(line, prefix)
                    if category == "info_lines" and item_id in canonical_info_ids:
                        continue
                    identity = item_id if item_id is not None else line
                    if identity in seen[category]:
                        continue
                    seen[category].add(identity)
                    collected[category].append(line)

        collected["safe_meta_lines"] = merge_safe_global_meta_lines(
            safe_meta_sources
        )
        return collected

    def _sample_id_from_data(self, sample_data, format_keys: List[str]) -> Optional[str]:
        """Extract the original source ID from one sample/evidence block.

        This is used to rebuild INFO/SOURCE_IDS in the same order as
        INFO/SOURCES and the retained sample columns.

        Raw FORMAT strings are parsed with the shared SVCF block parser so
        colon-containing source IDs and ALT values are preserved intact.
        """
        if sample_data is None:
            return None

        if isinstance(sample_data, dict):
            # Evidence-block ID is the source record ID for this block.
            # ``original_id`` is a compatibility fallback for older objects.
            source_id = sample_data.get("ID", sample_data.get("original_id"))
            if source_id not in (None, "", ".", "unknown"):
                return str(source_id)
            return None

        if "ID" not in format_keys:
            return None

        parsed = parse_svcf_sample_block(
            format_keys,
            str(sample_data),
        )
        value = parsed.get("ID")

        if value not in (None, "", ".", "unknown"):
            return str(value)

        return None

    def _write_event(self, file_handle, event):
        """Write one sample-mode SVCF record.

        event.ordered_samples must already be aligned to the global sample
        column order. This method does not infer sample identity; it only writes
        the fixed ordered columns.
        """
        all_sample_names = self.name_mapper.get_all_display_names()
        ordered_samples = getattr(event, "ordered_samples", [])

        # SVCF 1.1 sample mode has its own synthesized-call schema.
        # Keep ID:SC:REF:ALT:CO as the final five fields so the shared parser
        # remains robust to colon-containing IDs and BND/symbolic ALT values.
        format_field = SAMPLE_FORMAT_V11
        format_keys = format_field.split(":")

        sources = []
        source_ids = []

        for index, sample_data in enumerate(ordered_samples):
            if sample_data is not None:
                sources.append(all_sample_names[index])

                source_id = self._sample_id_from_data(sample_data, format_keys)
                source_ids.append(source_id if source_id else ".")

        # SVCF 1.1 positional INFO lists have no escaping layer.  Fail
        # loudly if a sample/source label or retained source ID contains an
        # INFO/list delimiter or whitespace rather than writing an ambiguous
        # file.
        sources = [validate_source_label(value) for value in sources]
        source_ids = [validate_source_id(value) for value in source_ids]

        sources_str = ",".join(sources) if sources else "."
        source_ids_str = ",".join(source_ids) if source_ids else "."

        # Prepare INFO field.
        # Do not inherit stale per-record SOURCES / SOURCE_IDS from the
        # representative event. They may already exist if the input is a
        # previously merged/subset/normalized SVCF. Always write exactly one
        # fresh SOURCES and one fresh SOURCE_IDS field for this output record.
        info_items = []
        for key, value in event.info.items():
            if key in {"SOURCES", "SOURCE_IDS"}:
                continue
            info_items.append(format_vcf_info_item(key, value))

        info_items.append(f"SOURCES={sources_str}")
        info_items.append(f"SOURCE_IDS={source_ids_str}")
        info_field = ";".join(info_items)

        sample_part = self._format_sample_columns(ordered_samples, format_keys)

        record = (
            f"{event.chrom}\t{event.pos}\t{event.sv_id}\t{event.ref}\t{event.alt}\t"
            f"{event.quality}\t{event.filter}\t{info_field}\t{format_field}\t"
            f"{sample_part}\n"
        )
        file_handle.write(record)

    def _format_sample_columns(self, ordered_samples, format_keys):
        """Format sample columns using preprocessed ordered_samples.

        Missing samples are represented as 0/0 placeholders with the same FORMAT
        field length.
        """
        sample_columns = []

        missing_values = {key: "." for key in format_keys}
        missing_values.update({"GT": "0/0", "AD": ".,.", "UC": "0", "UV": "0"})
        missing_data = ":".join(missing_values[key] for key in format_keys)

        for sample_data in ordered_samples:
            if sample_data is not None:
                if isinstance(sample_data, dict):
                    values = []
                    for key in format_keys:
                        value = sample_data.get(key, ".")
                        if isinstance(value, list):
                            value = ",".join(map(str, value))
                        elif value is None:
                            value = "."
                        values.append(str(value))
                    sample_str = ":".join(values)
                else:
                    sample_str = str(sample_data)

                sample_columns.append(sample_str)
            else:
                sample_columns.append(missing_data)

        return "\t".join(sample_columns)