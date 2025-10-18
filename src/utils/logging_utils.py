from __future__ import annotations

from datetime import datetime
import logging
import os
import pathlib
from typing import Any
from zoneinfo import ZoneInfo

from src.scripts.plaque_alignment.config import PathsCfg

# Global logger for this module
logger = logging.getLogger(__name__)


class TZFormatter(logging.Formatter):
    """
    Formatter that emits ISO-8601 timestamps in a specified IANA timezone
    (e.g., Europe/Zurich), including DST transitions.
    """

    def __init__(self, *args: Any, tz: str = "Europe/Zurich", **kwargs: Any):
        super().__init__(*args, **kwargs)
        self._tz = ZoneInfo(tz)

    def formatTime(self, record: logging.LogRecord, datefmt: str | None = None) -> str:  # noqa: N802
        # Use timezone-aware datetime derived from the record timestamp.
        # Example output: 2025-10-18T16:27:03.123+02:00
        dt = datetime.fromtimestamp(record.created, self._tz)
        if datefmt:
            # If a custom datefmt is supplied, honor it (still tz-aware)
            return dt.strftime(datefmt)
        # Default: ISO-8601 with milliseconds and numeric offset
        return dt.isoformat(timespec="milliseconds")


def _timestamp_for_filename(tzname: str = "Europe/Zurich") -> str:
    """Timestamp safe for filenames, e.g., 2025-10-18T16-27-03+02-00"""
    tz = ZoneInfo(tzname)
    # Use offset in filename so you can see whether it was CET (+01:00) or CEST (+02:00)
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

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(fmt)

    root.addHandler(fh)
    root.addHandler(ch)

    logger.debug("Logging configured (level=%s, tz=%s)", level_name, tzname)
    logger.info("Log file: %s", log_file)
    return log_file
