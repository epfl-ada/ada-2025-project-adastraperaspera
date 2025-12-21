from __future__ import annotations

from datetime import datetime
import logging
import os
import pathlib
import sys
from typing import Any
from zoneinfo import ZoneInfo

from src.scripts.plaque_alignment.config import PathsCfg


logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
    force=True,
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

if not logger.handlers:
    h = logging.StreamHandler(sys.stdout)
    h.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(h)

logger.propagate = False


class TZFormatter(logging.Formatter):
    """
    Formatter that emits ISO-8601 timestamps in a specified IANA timezone
    (e.g., Europe/Zurich), including DST transitions.
    """

    def __init__(self, *args: Any, tz: str = "Europe/Zurich", **kwargs: Any):
        super().__init__(*args, **kwargs)
        self._tz = ZoneInfo(tz)

    def formatTime(self, record: logging.LogRecord, datefmt: str | None = None) -> str:
        dt = datetime.fromtimestamp(record.created, self._tz)
        if datefmt:
            return dt.strftime(datefmt)
        return dt.isoformat(timespec="milliseconds")


class ImmediateFlushHandler(logging.StreamHandler):
    """StreamHandler that flushes immediately after each log record."""

    def __init__(self, stream=None):
        super().__init__(stream=sys.stdout if stream is None else stream)

    def emit(self, record):
        super().emit(record)
        self.flush()


class ImmediateFlushFileHandler(logging.FileHandler):
    """FileHandler that flushes immediately after each log record."""

    def emit(self, record):
        super().emit(record)
        self.flush()


def _timestamp_for_filename(tzname: str = "Europe/Zurich") -> str:
    """Timestamp safe for filenames, e.g., 2025-10-18T16-27-03+02-00"""
    tz = ZoneInfo(tzname)
    return (
        datetime.now(tz).strftime("%Y-%m-%dT%H-%M-%S%z")[:-2]
        + "-"
        + datetime.now(tz).strftime("%z")[-2:]
    )


def setup_logging(
    paths: PathsCfg, level_str: str | None = None, tzname: str = "Europe/Zurich"
) -> pathlib.Path:
    paths.logs_dir.mkdir(parents=True, exist_ok=True)
    run_ts = _timestamp_for_filename(tzname)
    log_file = paths.logs_dir / f"xenium_plaque_{run_ts}.log"

    level_name = (level_str or os.getenv("LOG_LEVEL") or "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)

    root = logging.getLogger()
    root.setLevel(level)
    for h in list(root.handlers):
        root.removeHandler(h)

    fmt = TZFormatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s", tz=tzname)

    fh = ImmediateFlushFileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    ch = ImmediateFlushHandler(stream=sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)

    root.addHandler(fh)
    root.addHandler(ch)

    logger.debug("Logging configured (level=%s, tz=%s)", level_name, tzname)
    root.info("Log file: %s", log_file)
    return log_file


def log_runtime(start_time: datetime, end_time: datetime) -> None:
    """Log the total runtime given start and end timestamps."""
    elapsed = end_time - start_time
    print(f"Total runtime: {elapsed:.1f} seconds")
