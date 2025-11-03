from __future__ import annotations

from pathlib import Path
import sys
import pytest

# Ensure `src` package is importable by adding project root to sys.path
# This supports imports like `from src.scripts...`.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def pytest_configure(config):
    """Configure pytest to suppress warnings from production code."""
    # Suppress UserWarning about non-interactive matplotlib backend (from plt.show() calls)
    config.addinivalue_line(
        "filterwarnings",
        "ignore:FigureCanvasAgg is non-interactive:UserWarning",
    )
    # Suppress DeprecationWarning about is_sparse (pandas deprecated API)
    config.addinivalue_line(
        "filterwarnings",
        "ignore:is_sparse is deprecated:DeprecationWarning",
    )
