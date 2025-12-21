from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("yaml")
from src.utils.logging_utils import setup_logging


def test_setup_logging_creates_log_file_and_writes(tmp_path: Path) -> None:
    logs_dir = tmp_path / "logs"
    paths_stub = SimpleNamespace(logs_dir=logs_dir)

    log_file = setup_logging(paths=paths_stub, level_str="DEBUG", tzname="UTC")

    assert log_file.parent == logs_dir
    assert log_file.exists(), "Log file should be created"

    logging.getLogger(__name__).info("hello from test")
    content = log_file.read_text(encoding="utf-8")
    assert "Log file:" in content
    assert "hello from test" in content
