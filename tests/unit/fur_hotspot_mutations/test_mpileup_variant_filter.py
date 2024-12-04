from collections import Counter
import logging
from pathlib import Path

import pandas as pd
import pytest

from fur_hotspot_mutations.mpileup_variant_filter import (
    tsv_to_df,
    extract_true_positives_from_mpileup_df,
    extract_false_negatives_from_mpileup_df,
    _check_tn_pairs,
    _extract_tn_pairs_from_df,
    _count_tn_pairs,
)
from tests.mocks.mock_files import (
    get_example_mpileup_file,
    get_example_tn_pairs_file,
)


# FIXTURES
@pytest.fixture
def mpileup_file() -> Path:
    """
    Fixture to provide an example mpileup file.

    Returns:
        Path: The path to the example mpileup file.
    """
    example_mpileup_file = get_example_mpileup_file()
    return example_mpileup_file


@pytest.fixture
def mpileup_df(mpileup_file) -> pd.DataFrame:
    """
    Fixture to provide a DataFrame of the example mpileup data.

    Returns:
        pd.DataFrame: The mpileup DataFrame.
    """
    df = pd.read_csv(mpileup_file, sep="\t", comment="#")
    return df


@pytest.fixture
def tn_pairs_file() -> Path:
    """
    Fixture to provide an example tumour-normal pairs file.

    Returns:
        Path: The path to the example tumour-normal pairs file.
    """
    example_tn_pairs_file = get_example_tn_pairs_file()
    return example_tn_pairs_file


@pytest.fixture
def tn_pairs_df(tn_pairs_file) -> pd.DataFrame:
    """
    Fixture to provide a DataFrame of the tumour-normal pairs data.

    Returns:
        pd.DataFrame: The tumour-normal pairs DataFrame.
    """
    df = pd.read_csv(tn_pairs_file, sep="\t")
    return df


# TESTS


# Tests for tsv_to_df


def test_tsv_to_df_valid_file(mpileup_file):
    """
    Test that tsv_to_df correctly reads a valid TSV file into a DataFrame.
    """
    # When
    df = tsv_to_df(mpileup_file)

    # Then
    assert isinstance(df, pd.DataFrame)
    assert not df.empty


def test_tsv_to_df_file_not_found(tmp_path):
    """
    Test that tsv_to_df raises a ValueError when the file does not exist.
    """
    # Given
    non_existent_file = tmp_path / "non_existent_file.tsv"

    # Then
    with pytest.raises(
        ValueError, match=r"Error: .* does not exist\. Please check the file path\."
    ):
        tsv_to_df(non_existent_file)


def test_tsv_to_df_empty_file(tmp_path):
    """
    Test that tsv_to_df raises a ValueError when the file is empty.
    """
    # Given
    empty_file = tmp_path / "empty_file.tsv"
    empty_file.write_text("")

    # Then
    with pytest.raises(ValueError, match=r"Error: .* is empty"):
        tsv_to_df(empty_file)


def test_tsv_to_df_invalid_format(tmp_path):
    """
    Test that tsv_to_df raises a ValueError when the file cannot be read due to invalid format.
    """
    # Given
    invalid_data = "This is not a valid TSV content"
    invalid_file = tmp_path / "invalid_file.tsv"
    invalid_file.write_text(invalid_data)

    # Then
    with pytest.raises(ValueError, match=r".* is not tab-delimited"):
        tsv_to_df(invalid_file)


# Tests for extract_true_positives_from_mpileup_df
def test_extract_true_positives_from_mpileup_df_valid_data(mpileup_df, caplog):
    """
    Test that extract_true_positives_from_mpileup_df correctly extracts true positives from valid data.
    """
    # Given
    caplog.set_level(logging.INFO)
    expected_true_pos_genes = ["TP53", "EGFR", "BRAF", "IDH1", "KRAS"]

    # When
    result_df = extract_true_positives_from_mpileup_df(mpileup_df)

    # Then
    assert list(result_df["Hugo_Symbol"]) == expected_true_pos_genes
    assert len(result_df) == len(expected_true_pos_genes)
    assert f"Detected {result_df.shape[0]} true positive variants." in caplog.text


def test_extract_true_positives_from_mpileup_df_no_true_pos(tmp_path, caplog):
    """
    Test that extract_true_positives_from_mpileup_df handles cases where there are no true positives.
    """
    # Given
    caplog.set_level(logging.INFO)
    no_true_pos_data = """Hugo_Symbol\tGene\tHGVSp_Short\tChromosome\tStart_Position\tEnd_Position\tVariant_Type\tReference_Allele\tTumour_Seq_Allele2\tTumor_Sample_Barcode\tAlt_Count\tTot_Count\tAlt_Perc\tVar_in_MAF\tStatus
TP53\tTP53\tp.R175H\t17\t7579472\t7579472\tSNP\tC\tT\tSample_T1\t45\t100\t0.45\tPRESENT_IN_MAF\tFALSE_NEGATIVE
"""
    no_true_pos_file = tmp_path / "no_true_pos.tsv"
    no_true_pos_file.write_text(no_true_pos_data)
    no_true_pos_df = pd.read_csv(no_true_pos_file, sep="\t", comment="#")

    # When
    result_df = extract_true_positives_from_mpileup_df(no_true_pos_df)

    # Then
    assert result_df.empty
    assert "No true positives detected" in caplog.text


# Tests for extract_false_negatives_from_mpileup_df
def test_extract_false_negatives_from_mpileup_df_valid_data(mpileup_df, caplog):
    """
    Test that extract_false_negatives_from_mpileup_df correctly extracts false negatives from valid data.
    """
    # Given
    caplog.set_level(logging.INFO)
    expected_false_neg_counts = Counter(
        {"TP53": 2, "EGFR": 3, "BRAF": 8, "IDH1": 3, "KRAS": 2}
    )

    # When
    result_df = extract_false_negatives_from_mpileup_df(mpileup_df)
    actual_false_neg_counts = Counter(result_df["Hugo_Symbol"])

    # Then
    assert actual_false_neg_counts == expected_false_neg_counts
    assert f"Detected {result_df.shape[0]} false negative variants." in caplog.text


def test_extract_false_negatives_from_mpileup_df_no_false_neg(tmp_path, caplog):
    """
    Test that extract_false_negatives_from_mpileup_df handles cases where there are no false negatives.
    """
    # Given
    caplog.set_level(logging.INFO)
    no_false_neg_data = """Hugo_Symbol\tGene\tHGVSp_Short\tChromosome\tStart_Position\tEnd_Position\tVariant_Type\tReference_Allele\tTumour_Seq_Allele2\tTumor_Sample_Barcode\tAlt_Count\tTot_Count\tAlt_Perc\tVar_in_MAF\tStatus
TP53\tTP53\tp.R175H\t17\t7579472\t7579472\tSNP\tC\tT\tSample_T1\t45\t100\t0.45\tPRESENT_IN_MAF\tTRUE_POSITIVE
"""
    no_false_neg_file = tmp_path / "no_false_neg.tsv"
    no_false_neg_file.write_text(no_false_neg_data)
    no_false_neg_df = pd.read_csv(no_false_neg_file, sep="\t", comment="#")

    # When
    result_df = extract_false_negatives_from_mpileup_df(no_false_neg_df)

    # Then
    assert result_df.empty
    assert "No false negatives detected" in caplog.text


# Tests for _check_tn_pairs
def test_check_tn_pairs_valid_data(mpileup_df, tn_pairs_df):
    """
    Test that _check_tn_pairs passes when mpileup and tumour-normal pairs have matching samples.
    """
    # Given
    # This test should pass as all samples in the example tumour-normal pairs are present in the mpileup DataFrame

    # When and Then
    _check_tn_pairs(mpileup_df, tn_pairs_df)


def test_check_tn_pairs_missing_sample_in_mpileup(mpileup_df, tn_pairs_df):
    """
    Test that _check_tn_pairs raises an error when a sample is missing in the mpileup DataFrame.
    """
    # Given
    sample_to_remove = "TUMOUR_SAMPLE_001"
    mpileup_df_filtered = mpileup_df[
        mpileup_df["Tumor_Sample_Barcode"] != sample_to_remove
    ]

    # When and Then
    with pytest.raises(
        ValueError, match=f"Missing samples in mpileup file: {sample_to_remove}"
    ):
        _check_tn_pairs(mpileup_df_filtered, tn_pairs_df)


def test_check_tn_pairs_missing_sample_in_tn_pairs(mpileup_df, tn_pairs_df):
    """
    Test that _check_tn_pairs raises an error when a sample is missing in the tumour-normal pairs DataFrame.
    """
    # Given
    sample_to_remove = "TUMOUR_SAMPLE_001"
    tn_pairs_df_filtered = tn_pairs_df[tn_pairs_df["TUMOUR"] != sample_to_remove]

    # When and Then
    with pytest.raises(
        ValueError,
        match=f"Missing samples in tumour-normal pair file: NORMAL_SAMPLE_001,TUMOUR_SAMPLE_001",
    ):
        _check_tn_pairs(mpileup_df, tn_pairs_df_filtered)


# Tests for _extract_tn_pairs_from_file
def test_extract_tn_pairs_from_df_valid_data(tn_pairs_df):
    """
    Test that _extract_tn_pairs_from_file correctly extracts tumour-normal pairs from valid data.
    """
    # Given
    expected_result = {
        "TUMOUR_SAMPLE_001_NORMAL_SAMPLE_001": {
            "TUMOUR": "TUMOUR_SAMPLE_001",
            "NORMAL": "NORMAL_SAMPLE_001",
        },
        "TUMOUR_SAMPLE_002_NORMAL_SAMPLE_002": {
            "TUMOUR": "TUMOUR_SAMPLE_002",
            "NORMAL": "NORMAL_SAMPLE_002",
        },
        "TUMOUR_SAMPLE_003_NORMAL_SAMPLE_003": {
            "TUMOUR": "TUMOUR_SAMPLE_003",
            "NORMAL": "NORMAL_SAMPLE_003",
        },
        # Add other expected pairs as needed
    }

    # When
    result = _extract_tn_pairs_from_df(tn_pairs_df)

    # Then
    for key in expected_result:
        assert key in result
        assert result[key] == expected_result[key]


# Tests for _count_tn_pairs
def test_count_tn_pairs(mpileup_df, tn_pairs_df):
    """
    Test that _count_tn_pairs correctly counts the number of tumour-normal pairs in the DataFrame.
    """
    # Given
    tn_pairs = _extract_tn_pairs_from_df(tn_pairs_df)
    # Assume we have a variant that exists in all tumour-normal pairs
    variant_rows_df = mpileup_df.copy()

    # When
    pair_count = _count_tn_pairs(variant_rows_df, tn_pairs)

    # Then
    assert pair_count == len(tn_pairs)


def test_count_tn_pairs_no_pairs(mpileup_df, tn_pairs_df):
    """
    Test that _count_tn_pairs returns zero when no pairs are present in the DataFrame.
    """
    # Given
    tn_pairs = _extract_tn_pairs_from_df(tn_pairs_df)
    # Create a DataFrame with no matching sample IDs
    variant_rows_df = mpileup_df[
        mpileup_df["Tumor_Sample_Barcode"] == "Nonexistent_Sample"
    ]

    # When
    pair_count = _count_tn_pairs(variant_rows_df, tn_pairs)

    # Then
    assert pair_count == 0
