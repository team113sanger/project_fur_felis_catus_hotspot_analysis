import subprocess
import shlex
import shutil

import fur_hotspot_mutations
import utils.constants

MODULE_NAME = fur_hotspot_mutations.__name__
PROGRAM_NAME = utils.constants.PROGRAM_NAME

# HELPERS


def get_subprocess_message(subproces_result: subprocess.CompletedProcess) -> str:
    indent = " " * 2

    msg = (
        f"Error running CLI command. "
        f"{indent}Command: {subproces_result.args}\n"
        f"{indent}Return code: {subproces_result.returncode}\n"
        f"{indent}Stdout: {subproces_result.stdout!r}\n"
        f"{indent}Stderr: {subproces_result.stderr!r}"
    )
    return msg


# TESTS


def test_python_dash_m__version():
    # Given
    cmd = f"python -m {MODULE_NAME} --version"
    expected_version = fur_hotspot_mutations.__version__

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout
    assert expected_version in subproces_result.stdout


def test_cli_on_path():
    # When
    result = shutil.which(PROGRAM_NAME)

    assert result is not None, f"{PROGRAM_NAME} is not in PATH, has the name changed?"


def test_cli__version():
    # Given
    cmd = f"{PROGRAM_NAME} --version"
    expected_version = fur_hotspot_mutations.__version__

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout
    assert expected_version in subproces_result.stdout


def test_cli__help():
    # Given
    cmd = f"{PROGRAM_NAME} --help"

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout


def test_cli__extract_hotspot_mutations__help():
    # Given
    command = "extract_hotspot_mutations"
    cmd = f"{PROGRAM_NAME} {command} --help"

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert command in subproces_result.stdout
