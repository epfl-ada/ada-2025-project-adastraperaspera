"""
Download and extract the Xenium TgCRND8 dataset.
"""

import os
import zipfile

import requests
from tqdm import tqdm

from src.utils.logging_utils import logger


def download_xenium_dataset(url: str, output_dir: str) -> str:
    """
    Download and extract the Xenium TgCRND8 dataset from 10x Genomics.

    Args:
        url (str): URL of the dataset ZIP file.
        output_dir (str): Directory where files should be extracted.

    Returns:
        str: Path to the directory containing the extracted files.
    """
    os.makedirs(output_dir, exist_ok=True)
    zip_path = os.path.join(output_dir, "Xenium_TgCRND8_17mo_outs.zip")

    logger.info(f"Downloading Xenium dataset from {url}")

    response = requests.get(url, stream=True, timeout=30)
    response.raise_for_status()
    total_size = int(response.headers.get("content-length", 0))

    with (
        open(zip_path, "wb") as file,
        tqdm(
            desc="Downloading Xenium dataset",
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar,
    ):
        for chunk in response.iter_content(1024):
            if chunk:
                size = file.write(chunk)
                bar.update(size)

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(output_dir)

    logger.info(f"Extracted dataset to {output_dir}")
    return output_dir
