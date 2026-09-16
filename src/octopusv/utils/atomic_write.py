"""Shared same-directory atomic output replacement.

Writers that may fail after output begins should write to the yielded temporary
path.  On successful context exit the temporary file replaces the target with
``os.replace``.  On any failure, including ``KeyboardInterrupt``, the temporary
file is removed and the previous target is left untouched.
"""

from __future__ import annotations

import os
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator


@contextmanager
def atomic_output_path(output_file: str | Path) -> Iterator[Path]:
    """Yield an exclusive same-directory temporary path and atomically commit it.

    Existing target permission bits are preserved.  New targets use normal
    text-file creation permissions (``0o666`` masked by the process umask).
    """
    output_path = Path(output_file)
    output_dir = output_path.parent

    existing_mode = None
    if output_path.exists():
        existing_mode = output_path.stat().st_mode & 0o777

    temp_path: Path | None = None
    fd: int | None = None

    try:
        for _ in range(10):
            candidate = output_dir / f".{output_path.name}.{uuid.uuid4().hex}.tmp"
            try:
                fd = os.open(
                    candidate,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                    0o666,
                )
                temp_path = candidate
                break
            except FileExistsError:
                continue

        if fd is None or temp_path is None:
            raise OSError(f"Could not create temporary output beside {output_path}")

        # Reserve the unique path, then let the caller reopen it using the API
        # its writer already expects (path-based or handle-based).  Preserve an
        # existing target mode only after writing: applying a read-only mode
        # before reopening the temporary path would make the writer fail.
        os.close(fd)
        fd = None

        yield temp_path

        if existing_mode is not None:
            os.chmod(temp_path, existing_mode)
        os.replace(temp_path, output_path)
        temp_path = None

    except BaseException:
        if fd is not None:
            os.close(fd)
        if temp_path is not None:
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass
        raise
