# This is where you may define fixtures that are shared across multiple test files.
#
# See this FAQ: https://docs.pytest.org/en/8.2.x/how-to/index.html
# Or read this section of the docs: https://docs.pytest.org/en/8.2.x/how-to/fixtures.html#scope-sharing-fixtures-across-classes-modules-packages-or-session

import os
from pathlib import Path
import importlib.resources

import pytest

import tests.mocks

# CONSTANTS
ENV_VAR_MAF_DIR = "TEST_MAF_DIR"


# FIXTURES
@pytest.fixture
def maf() -> Path:
    # If the environment variable is not set we use importlib.resources,
    # otherwise we use the path from the environment variable
    # to find the MAF file.
    file_name = "keepPA_vaf_size_filt_matched_6711_2820.maf"
    env_var_exists = bool(os.environ.get(ENV_VAR_MAF_DIR, None))

    if env_var_exists:
        from_env = True
        maf_file_dir = Path(os.environ[ENV_VAR_MAF_DIR])
        if not maf_file_dir.is_dir():
            msg = f"Environment variable {ENV_VAR_MAF_DIR!r} is not a directory: {maf_file_dir}"
            raise NotADirectoryError(msg)
        maf_file = maf_file_dir / file_name
    else:
        from_env = False
        resource = importlib.resources.files(tests.mocks) / file_name
        maf_file = Path(str(resource))  # Recast Resource obj to Path obj via str()

    # Check if the file exists
    if not maf_file.exists() and from_env:
        msg = f"Could not find MAF file at {maf_file} from environment variable {ENV_VAR_MAF_DIR!r}"
        raise FileNotFoundError(msg)
    if not maf_file.exists() and not from_env:
        msg = f"Could not find MAF file at {maf_file} via importlib.resources"
        raise FileNotFoundError(msg)

    return maf_file
