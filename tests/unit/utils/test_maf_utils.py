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


# Tests for is_maf_format
def test_is_maf_format(mock_file_path: Path):
    # Test with a valid MAF file
    valid_maf_data = (
        "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\t"
        "Variant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele2\n"
        "TP53\tTCGA\tGRCh38\t17\t7574003\t7574010\tMissense_Mutation\tSNP\tG\tA\n"
    )
    with patch("pathlib.Path.open", mock_open(read_data=valid_maf_data)):
        assert (
            is_maf_format(mock_file_path) is True
        ), "The function should return True for a valid MAF file."

    # Test with a file missing required columns
    missing_columns_data = "Some_Other_Header\tUnrelated_Column\n"
    with patch("pathlib.Path.open", mock_open(read_data=missing_columns_data)):
        assert (
            is_maf_format(mock_file_path) is False
        ), "The function should return False for a file missing required MAF columns."

    # Test with invalid Start_Position and End_Position values
    invalid_positions_data = (
        "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\t"
        "Variant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele2\n"
        "TP53\tTCGA\tGRCh38\t17\t7574010\t7574003\tMissense_Mutation\tSNP\tG\tA\n"
    )
    with patch("pathlib.Path.open", mock_open(read_data=invalid_positions_data)):
        assert (
            is_maf_format(mock_file_path) is False
        ), "The function should return False for invalid position data (Start_Position > End_Position)."

    # Test with non-integer Start_Position or End_Position
    non_integer_positions_data = (
        "Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\t"
        "Variant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele2\n"
        "TP53\tTCGA\tGRCh38\t17\tNaN\t7574010\tMissense_Mutation\tSNP\tG\tA\n"
    )
    with patch("pathlib.Path.open", mock_open(read_data=non_integer_positions_data)):
        assert (
            is_maf_format(mock_file_path) is False
        ), "The function should return False for non-integer positions."

    # Test with an empty file
    empty_file_data = ""
    with patch("pathlib.Path.open", mock_open(read_data=empty_file_data)):
        assert (
            is_maf_format(mock_file_path) is False
        ), "The function should return False for an empty file."

    # Test with a file containing only comments
    comment_only_data = "# This is a comment line\n# Another comment line\n"
    with patch("pathlib.Path.open", mock_open(read_data=comment_only_data)):
        assert (
            is_maf_format(mock_file_path) is False
        ), "The function should return False for a file with only comments."

    # Test with a file that raises an exception (e.g., unreadable file)
    with patch("pathlib.Path.open", side_effect=Exception("File read error")):
        assert (
            is_maf_format(mock_file_path) is False
        ), "The function should return False when an exception is raised while reading the file."

    # Test with a file missing header but having valid data
    no_header_data = (
        "TP53\tTCGA\tGRCh38\t17\t7574003\t7574010\tMissense_Mutation\tSNP\tG\tA\n"
    )
    with patch("pathlib.Path.open", mock_open(read_data=no_header_data)):
        assert (
            is_maf_format(mock_file_path) is False
        ), "The function should return False for a file with no valid header."


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


# Tests for check_column_consistency
def test_check_column_consistency(mock_file_path: Path):
    # Test with a valid file where all rows have consistent columns
    consistent_column_data = "col1\tcol2\tcol3\n1\t2\t3\n4\t5\t6\n7\t8\t9\n"
    with patch("pathlib.Path.open", mock_open(read_data=consistent_column_data)):
        assert (
            _check_column_consistency(mock_file_path) is True
        ), "The function should return True for a file with consistent column counts"

    # Test with a file with inconsistent columns
    inconsistent_column_data = "col1\tcol2\tcol3\n1\t2\t3\n4\t5\n7\t8\t9\n"
    with patch("pathlib.Path.open", mock_open(read_data=inconsistent_column_data)):
        assert (
            _check_column_consistency(mock_file_path) is False
        ), "The function should return False for a file with inconsistent column counts"

    # Test with an empty file
    empty_data = ""
    with patch("pathlib.Path.open", mock_open(read_data=empty_data)):
        assert (
            _check_column_consistency(mock_file_path) is False
        ), "The function should return False for an empty file"

    # Test with a file that does not exist
    non_existent_path = Path("nonexistent_file.txt")
    assert _check_column_consistency(non_existent_path) is False


def test_check_column_consistency_with_valid_file(maf):
    # Ensure the mock MAF file exists
    assert maf.exists(), "The mock MAF file does not exist."

    # Test the function with a valid MAF file
    assert (
        _check_column_consistency(maf) is True
    ), "The function should return True for a valid MAF file."
