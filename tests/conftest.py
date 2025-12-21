from __future__ import annotations

from pathlib import Path
import sys
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def pytest_configure(config):
    """Configure pytest to suppress warnings from production code."""

    config.addinivalue_line(
        "filterwarnings",
        "ignore:FigureCanvasAgg is non-interactive:UserWarning",
    )

    config.addinivalue_line(
        "filterwarnings",
        "ignore:is_sparse is deprecated:DeprecationWarning",
    )
