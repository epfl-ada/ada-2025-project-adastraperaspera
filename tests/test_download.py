"""
Unit tests for src.data.download module.
"""

from __future__ import annotations

import io
from pathlib import Path
import zipfile

import pytest

from src.data.download import download_xenium_dataset


class DummyResponse:
    """Lightweight stand-in for requests.get()"""

    def __init__(self, content: bytes):
        self.content = content
        self.headers = {"content-length": str(len(content))}

    def iter_content(self, chunk_size: int = 1024):
        yield self.content

    def raise_for_status(self):
        return None


@pytest.fixture
def fake_zip_bytes() -> bytes:
    """Create an in-memory ZIP file containing one dummy text file."""
    mem_zip = io.BytesIO()
    with zipfile.ZipFile(mem_zip, "w") as zf:
        zf.writestr("dummy.txt", "Hello ADA project!")
    mem_zip.seek(0)
    return mem_zip.read()


def test_download_and_extract(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fake_zip_bytes: bytes
) -> None:
    """Ensure the dataset is downloaded and extracted correctly (mocked)."""

    monkeypatch.setattr("requests.get", lambda *a, **kw: DummyResponse(fake_zip_bytes))

    out_dir = download_xenium_dataset("http://fake.url/data.zip", str(tmp_path))

    out_dir_path = Path(out_dir)
    extracted_files = list(out_dir_path.glob("*"))
    assert any(f.name == "dummy.txt" for f in extracted_files)
    assert out_dir_path.exists()
    assert out_dir_path.is_dir()

    with open(out_dir_path / "dummy.txt") as f:
        assert "Hello ADA project!" in f.read()


def test_download_raises_on_bad_url(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Ensure that failed requests raise appropriate errors."""

    class FailingResponse:
        def raise_for_status(self):
            raise Exception("HTTP error!")

    monkeypatch.setattr("requests.get", lambda *a, **kw: FailingResponse())

    with pytest.raises(Exception, match="HTTP error"):
        download_xenium_dataset("http://bad.url", str(tmp_path))
