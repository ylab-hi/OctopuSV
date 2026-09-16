"""Small shared helpers for transparent text input.

Compression is detected from the gzip magic bytes rather than the filename, so
plain text named ``*.gz`` and gzip/bgzip input without a compressed suffix are
both handled correctly.  Output compression is intentionally outside this
helper.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import IO


_GZIP_MAGIC = b"\x1f\x8b"


def is_gzip_path(path: str | Path) -> bool:
    """Return whether *path* contains gzip-compatible data.

    The historical function name is retained for callers, but detection is by
    file content, not suffix.
    """
    # Use os.open/os.read rather than the text reader itself so sniffing stays
    # independent of any wrapped/streaming builtins.open implementation.
    fd = os.open(os.fspath(path), os.O_RDONLY)
    try:
        return os.read(fd, 2) == _GZIP_MAGIC
    finally:
        os.close(fd)


def open_text_auto(
    path: str | Path,
    mode: str = "rt",
    *,
    encoding: str = "utf-8-sig",
) -> IO[str]:
    """Open plain or gzip/bgzip text input transparently.

    ``gzip.open`` supports concatenated gzip members, which also covers normal
    sequential reading of bgzip files without requiring a bgzip executable.
    """
    if "b" in mode:
        raise ValueError("open_text_auto only supports text modes.")
    if "r" not in mode:
        raise ValueError("open_text_auto is for input/read modes only.")

    text_mode = mode if "t" in mode else mode + "t"
    if is_gzip_path(path):
        return gzip.open(path, text_mode, encoding=encoding)
    return open(path, text_mode, encoding=encoding)
