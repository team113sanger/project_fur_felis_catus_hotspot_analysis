import pytest
import re

# https://gitlab.internal.sanger.ac.uk/help/user/packages/pypi_repository/index.md#ensure-your-version-string-is-valid
SEMVER_REGEX = "".join(
    [
        r"(?:",
        r"(?:([0-9]+)!)?",  # epoch
        r"([0-9]+(?:\.[0-9]+)*)",  # release segment
        r"([-_\.]?((a|b|c|rc|alpha|beta|pre|preview))[-_\.]?([0-9]+)?)?",  # pre-release
        r"((?:-([0-9]+))|(?:[-_\.]?(post|rev|r)[-_\.]?([0-9]+)?))?",  # post release
        r"([-_\.]?(dev)[-_\.]?([0-9]+)?)?",  # dev release
        r"(?:\+([a-z0-9]+(?:[-_\.][a-z0-9]+)*))?",  # local version
        r")",
    ]
)


@pytest.mark.filterwarnings("ignore")
def test_package_structure_by_importing_from_src():
    from src import fur_hotspot_mutations as module

    assert bool(module)


def test_package_structure_by_importing_from_installed_pkg():
    try:
        import fur_hotspot_mutations as module
    except ImportError:
        try:
            from src import fur_hotspot_mutations as module
        except ImportError as err:
            msg = "Could not import entrypoints from either src/ or from installed package. Something bad has happened."
            raise RuntimeError(msg) from err
        else:
            msg = "Could not import from package, but could import from src/\nDid you forget to install the package?\nEither do `poetry install` or `pip install -e .`"
    else:
        msg = ""
    assert bool(module), msg


def test_version():
    try:
        import fur_hotspot_mutations as module
    except ImportError:
        should_raise = True
        version = None
    else:
        should_raise = False
        version = module.__version__

    if should_raise:
        assert False, "Could not import package. Can't test version."

    err_msg = "Version is 0.0.0. Did you forget to set the version in pyproject.toml"
    assert version != "0.0.0" and version is not None, err_msg


def test_version_is_compatible():
    import fur_hotspot_mutations as module

    # Given
    rgx_pattern = re.compile(SEMVER_REGEX, re.VERBOSE)

    # When
    version = module.__version__

    # Then
    msg = f"Version {version} is not compatible with PEP440"
    assert rgx_pattern.match(version), msg
