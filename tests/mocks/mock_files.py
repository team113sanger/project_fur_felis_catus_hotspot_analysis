from importlib import resources
from pathlib import Path


def get_example_variant_file() -> Path:
    file = resources.files("tests.mocks").joinpath("example_variant_file.tsv")
    as_path = Path(str(file))
    if not as_path.exists():
        raise FileNotFoundError(f"File not found: {file}")

    return as_path


def get_example_maf() -> Path:
    file = resources.files("tests.mocks").joinpath("example.maf")
    as_path = Path(str(file))
    if not as_path.exists():
        raise FileNotFoundError(f"File not found: {file}")

    return as_path


def get_example_mpileup_file() -> Path:
    file = resources.files("tests.mocks").joinpath("example_mpileup_file.tsv")
    as_path = Path(str(file))
    if not as_path.exists():
        raise FileNotFoundError(f"File not found: {file}")

    return as_path


def get_example_tn_pairs_file() -> Path:
    file = resources.files("tests.mocks").joinpath("example_tn_pairs.tsv")
    as_path = Path(str(file))
    if not as_path.exists():
        raise FileNotFoundError(f"File not found: {file}")

    return as_path
