import os
from typing import List, Optional


class NameMapper:
    """Handle file-to-display-name mapping for both caller and sample modes."""

    def __init__(self, input_files: List[str], mode: str = "caller", custom_names: Optional[List[str]] = None):
        """
        Initialize name mapper.

        Args:
            input_files: List of input file paths
            mode: Either "caller" or "sample"
            custom_names: Optional custom names provided by user
        """
        self.input_files = [str(f) for f in input_files]
        self.mode = mode
        self.custom_names = custom_names
        self._validate_custom_names()

    def _validate_custom_names(self):
        """Validate that custom names count matches file count."""
        if self.custom_names and len(self.custom_names) != len(self.input_files):
            raise ValueError(
                f"Number of custom names ({len(self.custom_names)}) doesn't match number of files ({len(self.input_files)})")

    def get_display_name(self, file_path: str) -> str:
        """Get display name for a single file."""
        if self.custom_names:
            try:
                file_index = self.input_files.index(str(file_path))
                return self.custom_names[file_index]
            except (ValueError, IndexError):
                pass

        # Default: extract from filename
        return os.path.splitext(os.path.basename(str(file_path)))[0]

    def get_all_display_names(self) -> List[str]:
        """Get display names for all input files."""
        return [self.get_display_name(f) for f in self.input_files]

    def convert_source_string(self, source_file_str: str) -> str:
        """Convert comma-separated source file string to display names."""
        source_files = source_file_str.split(",")
        display_names = []
        for source_file in source_files:
            display_names.append(self.get_display_name(source_file.strip()))
        return ",".join(display_names)