import os


def normalize_source_path(source_file):
    """Return a stable real-path identity for one source/input file."""
    return os.path.normcase(
        os.path.realpath(
            os.path.abspath(str(source_file))
        )
    )
