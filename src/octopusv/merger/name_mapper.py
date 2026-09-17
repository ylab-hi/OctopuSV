import os
from typing import List, Optional


_COMPRESSED_SUFFIXES = (".gz", ".bgz", ".bgzip")


def default_label_from_path(file_path: str) -> str:
    """Derive a stable default display label from one input path.

    Compression should not change source/sample identity, so ``sample.svcf``
    and ``sample.svcf.gz`` both map to ``sample``.
    """
    name = os.path.basename(str(file_path))
    lowered = name.lower()
    for suffix in _COMPRESSED_SUFFIXES:
        if lowered.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return os.path.splitext(name)[0]


class NameMapper:
    """Handle file-to-display-name mapping for both caller and sample modes."""

    def __init__(self, input_files: List[str], mode: str = "caller", custom_names: Optional[List[str]] = None):
        self.input_files = [str(f) for f in input_files]
        self.mode = mode
        self.custom_names = custom_names
        self._validate_custom_names()

    def _validate_custom_names(self):
        if self.custom_names and len(self.custom_names) != len(self.input_files):
            raise ValueError(
                f"Number of custom names ({len(self.custom_names)}) doesn't match number of files ({len(self.input_files)})"
            )

    def get_display_name(self, file_path: str) -> str:
        if self.custom_names:
            try:
                file_index = self.input_files.index(str(file_path))
            except ValueError as exc:
                raise ValueError(
                    f"Input path {str(file_path)!r} is not present in the NameMapper input list; "
                    "refusing to fall back to a filename-derived label while custom names are active."
                ) from exc
            return self.custom_names[file_index]
        return default_label_from_path(file_path)

    def get_all_display_names(self) -> List[str]:
        return [self.get_display_name(f) for f in self.input_files]

    def convert_source_string(self, source_file_str: str) -> str:
        source_files = source_file_str.split(",")
        display_names = []
        for source_file in source_files:
            display_names.append(self.get_display_name(source_file.strip()))
        return ",".join(display_names)
