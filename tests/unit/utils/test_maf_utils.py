import typing as t
import logging
import os
from pathlib import Path
from unittest.mock import mock_open, patch

import pytest

from utils.maf_utils import (
    is_maf_format,
    _check_for_duplicates,
    _check_column_consistency,
)

from utils import constants

# TESTS


@pytest.mark.parametrize(
    "file_data, expected_result",
    [
        pytest.param(
            (
                "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\t"
                "Variant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele2\n"
                "TP53\tTCGA\tGRCh38\t17\t7574003\t7574010\tMissense_Mutation\tSNP\tG\tA\n"
            ),
            True,
            id="valid_maf_file",
        ),
        pytest.param(
            "Some_Other_Header\tUnrelated_Column\n",
            False,
            id="missing_required_columns",
        ),
        pytest.param(
            (
                "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\t"
                "Variant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele2\n"
                "TP53\tTCGA\tGRCh38\t17\t7574010\t7574003\tMissense_Mutation\tSNP\tG\tA\n"
            ),
            False,
            id="invalid_positions",
        ),
        pytest.param(
            (
                "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\t"
                "Variant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele2\n"
                "TP53\tTCGA\tGRCh38\t17\tNaN\t7574010\tMissense_Mutation\tSNP\tG\tA\n"
            ),
            False,
            id="non_integer_positions",
        ),
        pytest.param(
            "",
            False,
            id="empty_file",
        ),
        pytest.param(
            "# This is a comment line\n# Another comment line\n",
            False,
            id="comments_only",
        ),
        pytest.param(
            "TP53\tTCGA\tGRCh38\t17\t7574003\t7574010\tMissense_Mutation\tSNP\tG\tA\n",
            False,
            id="no_header",
        ),
    ],
)
def test_is_maf_format(
    tmp_path: Path,
    file_data: str,
    expected_result: bool,
):
    # Given
    mock_maf_file = tmp_path / "mock_maf_file.txt"
    mock_maf_file.write_text(file_data)
    assert mock_maf_file.exists(), "Precondition"

    # When
    result = is_maf_format(mock_maf_file)

    # Then
    assert result == expected_result


def test_is_maf_format__with_corrupted_file(tmp_path: Path):
    # Test with a file that raises an exception (e.g., unreadable file)
    mock_maf_file = tmp_path / "mock_maf_file.txt"
    with patch("pathlib.Path.open", side_effect=Exception("File read error")):
        assert (
            is_maf_format(mock_maf_file) is False
        ), "The function should return False when an exception is raised while reading the file."


def test_is_maf_format_with_valid_file(maf):
    # Ensure the mock MAF file exists
    assert maf.exists(), "The mock MAF file does not exist."

    # Test the function with a valid MAF file
    assert (
        is_maf_format(maf) is True
    ), "The function should return True for a valid MAF file."


# Tests for _check_for_duplicates
@pytest.mark.parametrize(
    "file_data, expected_result, expected_logging_warnings",
    [
        pytest.param(
            "line1\nline2\nline3\n# comment line\nline4\n",
            True,
            [],
            id="no_duplicates",
        ),
        pytest.param(
            "line1\nline2\nline3\nline1\nline4\n",
            False,
            ["Duplicate line found", "Line 1 and Line 4"],
            id="has_duplicates",
        ),
        pytest.param(
            "# Comment 1\n# Comment 2\n",
            False,
            ["contains no valid lines"],
            id="no_duplicates#comments_only",
        ),
        pytest.param(
            "",
            False,
            ["contains no valid lines"],
            id="empty_file",
        ),
    ],
)
def test_check_for_duplicates(
    tmp_path: Path,
    file_data: str,
    expected_result: bool,
    expected_logging_warnings: list[str],
    caplog,
):
    # Given
    mock_maf_file = tmp_path / "mock_maf_file.txt"
    mock_maf_file.write_text(file_data)
    assert mock_maf_file.exists(), "Precondition"

    # When
    with caplog.at_level(logging.WARNING, logger=constants.LOGGER_NAME):
        result = _check_for_duplicates(mock_maf_file)

    # Then
    assert result == expected_result
    for warning in expected_logging_warnings:
        assert warning in caplog.text, f"Expected warning '{warning}' not found in log."


def test_check_for_duplicates_with_valid_file(maf):
    # Ensure the mock MAF file exists
    assert maf.exists(), "The mock MAF file does not exist."

    # Test the function with a valid MAF file
    assert (
        _check_for_duplicates(maf) is True
    ), "The function should return True for a valid MAF file."


@pytest.mark.parametrize(
    "file_data, expected_result",
    [
        pytest.param(
            "col1\tcol2\tcol3\n1\t2\t3\n4\t5\t6\n7\t8\t9\n",
            True,
            id="consistent_columns",
        ),
        pytest.param(
            "col1\tcol2\tcol3\n1\t2\t3\n4\t5\n7\t8\t9\n",
            False,
            id="inconsistent_columns",
        ),
        pytest.param(
            "",
            False,
            id="empty_file",
        ),
    ],
)
def test_check_column_consistency(
    tmp_path: Path,
    file_data: str,
    expected_result: bool,
):
    # Given
    mock_file = tmp_path / "mock_column_file.txt"
    mock_file.write_text(file_data or "")
    assert mock_file.exists(), "Precondition"

    # When
    result = _check_column_consistency(mock_file)

    # Then
    assert result == expected_result


def test_check_column_consistency__non_existent_file():
    # Given
    non_existent_file = Path("non_existent_file.txt")
    assert not non_existent_file.exists(), "Precondition: The file should not exist."

    # When
    result = _check_column_consistency(non_existent_file)

    # Then
    assert result is False, "The function should return False for a non-existent file."


def test_check_column_consistency__with_valid_file(maf: Path):
    # Ensure the mock MAF file exists
    assert maf.exists(), "The mock MAF file does not exist."

    # Test the function with a valid MAF file
    assert (
        _check_column_consistency(maf) is True
    ), "The function should return True for a valid MAF file."
