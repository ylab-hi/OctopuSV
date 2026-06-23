import datetime
from typing import List, Optional

from .name_mapper import NameMapper


class MultiSampleWriter:
    """Write SVCF output for sample mode with fixed sample-column ordering.

    In sample mode, the #CHROM trailing columns are true sample/input columns.
    Their order is fixed by the input file order or by --sample-names.

    This writer assumes each event has event.ordered_samples already prepared
    in that exact order.
    """

    def __init__(self, name_mapper: NameMapper):
        """Initialize writer with a NameMapper."""
        self.name_mapper = name_mapper

    def write_results(self, output_file, events, contigs):
        """Write merged sample-mode SVCF results.

        Args:
            output_file: Output SVCF path.
            events: Merged events with ordered_samples attached.
            contigs: Dictionary of contig ID -> length.
        """
        with open(output_file, "w") as f:
            self._write_header(f, contigs)

            for event in events:
                self._write_event(f, event)

    def _write_header(self, file_handle, contigs):
        """Write SVCF header for sample mode."""
        file_handle.write("##fileformat=VCFv4.2\n")
        file_handle.write("##OctopuSV_mode=multi\n")

        file_date = datetime.datetime.now().strftime("%Y-%m-%d|%I:%M:%S%p|")
        file_handle.write(f"##fileDate={file_date}\n")
        file_handle.write("##source=OctopuSV\n")
        file_handle.write(
            "##OctopuSV_WARNING=This is SVCF format. "
            "Use 'octopusv svcf2vcf' to change back to standard VCF format "
            "before bcftools/vcftools\n"
        )

        for contig_id, contig_length in contigs.items():
            file_handle.write(f"##contig=<ID={contig_id},length={contig_length}>\n")

        self._write_standard_definitions(file_handle)

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
            '##INFO=<ID=AF,Number=1,Type=String,Description="Allele frequency">\n'
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
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Length of SV">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand orientation of SV">\n'
        )
        file_handle.write(
            '##FORMAT=<ID=QV,Number=1,Type=Integer,Description="Quality value">\n'
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

    def _sample_id_from_data(self, sample_data, format_keys: List[str]) -> Optional[str]:
        """Extract the original source ID from one sample/evidence block.

        This is used to rebuild INFO/SOURCE_IDS in the same order as
        INFO/SOURCES and the retained sample columns.

        For raw FORMAT strings, use a limited left split instead of splitting the
        whole string. This avoids breaking on ':' inside BND ALT values.
        """
        if sample_data is None:
            return None

        if isinstance(sample_data, dict):
            original_id = sample_data.get("original_id", sample_data.get("ID"))
            if original_id not in (None, "", ".", "unknown"):
                return str(original_id)
            return None

        if "ID" not in format_keys:
            return None

        id_index = format_keys.index("ID")

        # Split only enough fields to reach ID. FORMAT fields before ID should
        # not contain ':' in OctopuSV SVCF, while later ALT/CO fields may.
        values = str(sample_data).split(":", id_index + 1)

        if id_index >= len(values):
            return None

        value = values[id_index]
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

        format_field = event.format
        format_keys = format_field.split(":")

        sources = []
        source_ids = []

        for index, sample_data in enumerate(ordered_samples):
            if sample_data is not None:
                sources.append(all_sample_names[index])

                source_id = self._sample_id_from_data(sample_data, format_keys)
                source_ids.append(source_id if source_id else ".")

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
            info_items.append(f"{key}={value}")

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

        missing_data = "0/0" + ":" + ":".join(["."] * (len(format_keys) - 1))

        for sample_data in ordered_samples:
            if sample_data is not None:
                if isinstance(sample_data, dict):
                    values = [str(sample_data.get(key, ".")) for key in format_keys]
                    sample_str = ":".join(values)
                else:
                    sample_str = str(sample_data)

                if sample_str.endswith(":.:."):
                    sample_str = sample_str[:-4]

                sample_columns.append(sample_str)
            else:
                sample_columns.append(missing_data)

        return "\t".join(sample_columns)