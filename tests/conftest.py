# This is where you may define fixtures that are shared across multiple test files.
#
# See this FAQ: https://docs.pytest.org/en/8.2.x/how-to/index.html
# Or read this section of the docs: https://docs.pytest.org/en/8.2.x/how-to/fixtures.html#scope-sharing-fixtures-across-classes-modules-packages-or-session

import os
from pathlib import Path

import pytest

# CONSTANTS
ENV_VAR_MAF_DIR = "TEST_MAF_DIR"


# FIXTURES
@pytest.fixture
def mock_file_path(tmp_path):
    """Fixture to create a temporary file for testing."""
    return tmp_path / "test_file"


@pytest.fixture
def maf() -> Path:
    raw_maf_dir = os.environ.get(ENV_VAR_MAF_DIR, "tests/mocks")
    maf_dir = Path(raw_maf_dir)
    expected_maf = "keepPA_vaf_size_filt_matched_6711_2820.maf"
    maf = maf_dir / expected_maf
    if not maf.exists():
        raise FileNotFoundError(f"Could not find MAF file at {str(maf)}")

    return maf
