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
    _count_germline_tn_pairs,
    _get_matched_normal,
    _convert_chromosome_to_string,
    _construct_variant_file_row,
    write_variants_to_variant_file,
    process_true_positives,
    process_false_negatives,
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
    }

    # When
    result = _extract_tn_pairs_from_df(tn_pairs_df)

    # Then
    for key in expected_result:
        assert key in result
        assert result[key] == expected_result[key]


# # Tests for _count_tn_pairs
# def test_count_tn_pairs(mpileup_df, tn_pairs_df):
#     """
#     Test that _count_tn_pairs correctly counts the number of tumour-normal pairs in the DataFrame.
#     """
#     # Given
#     tn_pairs = _extract_tn_pairs_from_df(tn_pairs_df)
#     # Assume we have a variant that exists in all tumour-normal pairs
#     variant_rows_df = mpileup_df.copy()

#     # When
#     pair_count = _count_tn_pairs(variant_rows_df, tn_pairs)

#     # Then
#     assert pair_count == len(tn_pairs)


# def test_count_tn_pairs_no_pairs(mpileup_df, tn_pairs_df):
#     """
#     Test that _count_tn_pairs returns zero when no pairs are present in the DataFrame.
#     """
#     # Given
#     tn_pairs = _extract_tn_pairs_from_df(tn_pairs_df)
#     # Create a DataFrame with no matching sample IDs
#     variant_rows_df = mpileup_df[
#         mpileup_df["Tumor_Sample_Barcode"] == "Nonexistent_Sample"
#     ]

#     # When
#     pair_count = _count_tn_pairs(variant_rows_df, tn_pairs)

#     # Then
#     assert pair_count == 0


# # Tests for _count_alt_reads_in_normals


# def test_count_alt_reads_in_normals():
#     # Sample variant rows DataFrame
#     data = {
#         "Hugo_Symbol": ["TP53", "TP53", "TP53", "TP53", "TP53", "TP53"],
#         "Chromosome": ["chr17", "chr17", "chr17", "chr17", "chr17", "chr17"],
#         "Start_Position": [7579472, 7579472, 7579472, 7579472, 7579472, 7579472],
#         "Tumor_Sample_Barcode": [
#             "TUMOUR_SAMPLE_1",
#             "NORMAL_SAMPLE_1",
#             "TUMOUR_SAMPLE_2",
#             "NORMAL_SAMPLE_2",
#             "TUMOUR_SAMPLE_3",
#             "NORMAL_SAMPLE_3",
#         ],
#         "Alt_Count": [10, 5, 8, 2, 15, 4],
#         "Status": [
#             "FALSE_NEGATIVE",
#             "TRUE_NEGATIVE",
#             "FALSE_NEGATIVE",
#             "TRUE_NEGATIVE",
#             "FALSE_NEGATIVE",
#             "TRUE_NEGATIVE",
#         ],
#     }
#     variant_rows_df = pd.DataFrame(data)

#     # Sample tumor-normal pairs
#     tn_pairs = {
#         "PAIR_1": {"TUMOUR": "TUMOUR_SAMPLE_1", "NORMAL": "NORMAL_SAMPLE_1"},
#         "PAIR_2": {"TUMOUR": "TUMOUR_SAMPLE_2", "NORMAL": "NORMAL_SAMPLE_2"},
#         "PAIR_3": {"TUMOUR": "TUMOUR_SAMPLE_3", "NORMAL": "NORMAL_SAMPLE_3"},
#     }

#     # Set the minimum ALT reads threshold
#     min_alt_norm_reads = 3

#     # Expected count:
#     # NORMAL_SAMPLE_1 has Alt_Count = 5 (>3) --> counts
#     # NORMAL_SAMPLE_2 has Alt_Count = 2 (<=3) --> does not count
#     # NORMAL_SAMPLE_3 has Alt_Count = 4 (>3) --> counts
#     # Expected count = 2

#     # Call the function
#     count = _count_alt_reads_in_normals(variant_rows_df, tn_pairs, min_alt_norm_reads)

#     # Assert the expected count
#     assert count == 2, f"Expected count of 2, but got {count}"


# Tests for _count_germline_tn_pairs


def test_count_germline_tn_pairs():
    # Given
    variant_row_data = {
        "Tumor_Sample_Barcode": ["TUM1", "NORM1", "TUM2", "NORM2"],
        "Status": [
            "TRUE_NEGATIVE",
            "FALSE_NEGATIVE",
            "FALSE_NEGATIVE",
            "FALSE_NEGATIVE",
        ],
        "Alt_Count": [0, 10, 3, 12],  # Simulated ALT reads for tumour and normal
    }
    variant_row_df = pd.DataFrame(variant_row_data)

    tn_pairs = {
        "PAIR1": {"TUMOUR": "TUM1", "NORMAL": "NORM1"},
        "PAIR2": {"TUMOUR": "TUM2", "NORMAL": "NORM2"},
    }
    min_alt_norm_reads = 8

    # When
    germline_pair_count = _count_germline_tn_pairs(
        variant_row_df, tn_pairs, min_alt_norm_reads
    )

    # Then
    assert (
        germline_pair_count == 2
    ), f"Expected 2 germline pairs, got {germline_pair_count}"


# Tests for _get_matched_normal


def test_get_matched_normal_found(tn_pairs_df):
    tn_pairs = _extract_tn_pairs_from_df(tn_pairs_df)

    tumor_sample = "TUMOUR_SAMPLE_001"
    expected_normal = "NORMAL_SAMPLE_001"
    assert _get_matched_normal(tumor_sample, tn_pairs) == expected_normal


def test_get_matched_normal_not_found(tn_pairs_df):
    tn_pairs = _extract_tn_pairs_from_df(tn_pairs_df)

    tumor_sample = "NotASample"  # This tumor sample does not exist in the pairs
    with pytest.raises(ValueError, match="Error when extracting matched normal"):
        _get_matched_normal(tumor_sample, tn_pairs)


def test_get_matched_normal_empty_dict():
    tumor_sample = "TumorA"
    tn_pairs = {}  # Empty dictionary
    with pytest.raises(ValueError, match="Error when extracting matched normal"):
        _get_matched_normal(tumor_sample, tn_pairs)


# Tests for _convert_chromosome_to_string
def test_convert_chromosome_to_string():
    # Test with integer input
    assert _convert_chromosome_to_string(1) == "1"
    assert _convert_chromosome_to_string(23) == "23"

    # Test with string input
    assert _convert_chromosome_to_string("X") == "X"
    assert _convert_chromosome_to_string("Y") == "Y"

    # Test with edge cases
    with pytest.raises(ValueError, match="Value for chromosome is empty."):
        _convert_chromosome_to_string("")


# Tests for _construct_variant_file_row
def test_construct_variant_file_row(tn_pairs_df):
    # Given
    variant_row = pd.Series(
        {
            "Tumor_Sample_Barcode": "TUMOUR_SAMPLE_001",
            "Hugo_Symbol": "TP53",
            "Chromosome": 17,
            "Start_Position": 100000,
            "End_Position": 200000,
            "Reference_Allele": "C",
            "Tumour_Seq_Allele2": "T",
        }
    )
    tn_pairs = {
        "TUMOUR_SAMPLE_001_tNORMAL_SAMPLE_001": {
            "TUMOUR": "TUMOUR_SAMPLE_001",
            "NORMAL": "NORMAL_SAMPLE_001",
        }
    }
    add_to_maf = True

    expected_row = (
        "TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001\tTP53\t17\t100000\t200000\tC\tT\tADD\n"
    )

    # When
    result_row = _construct_variant_file_row(variant_row, tn_pairs, add_to_maf)

    # Then
    assert result_row == expected_row


# Test with missing field
def test_construct_variant_file_row_missing_field(tn_pairs_df):
    # Given
    variant_row = pd.Series(
        {
            "Tumor_Sample_Barcode": "TUMOUR_SAMPLE_001",
            "Hugo_Symbol": "TP53",
            "Chromosome": 17,
            "Start_Position": 100000,
            "End_Position": 200000,
            "Reference_Allele": "C"
            #'Tumour_Seq_Allele2': 'T'
        }
    )
    tn_pairs = {
        "TUMOUR_SAMPLE_001_tNORMAL_SAMPLE_001": {
            "TUMOUR": "TUMOUR_SAMPLE_001",
            "NORMAL": "NORMAL_SAMPLE_001",
        }
    }
    add_to_maf = True

    # When and Then
    with pytest.raises(ValueError, match="Variant row is missing one or more fields"):
        _construct_variant_file_row(variant_row, tn_pairs, add_to_maf)


def test_write_variants_to_variant_file(tmp_path):
    # Given
    variant_file = tmp_path / "variants_output.tsv"

    # A single variant row with required fields
    variant_row = pd.Series(
        {
            "Tumor_Sample_Barcode": "TUMOUR_SAMPLE_001",
            "Hugo_Symbol": "TP53",
            "Chromosome": "17",
            "Start_Position": 7579472,
            "End_Position": 7579472,
            "Reference_Allele": "C",
            "Tumour_Seq_Allele2": "T",
        }
    )

    # Tumour-normal pairs dictionary, must include a pair for "TUMOUR_SAMPLE_001"
    tn_pairs = {
        "TUMOUR_SAMPLE_001_NORMAL_SAMPLE_001": {
            "TUMOUR": "TUMOUR_SAMPLE_001",
            "NORMAL": "NORMAL_SAMPLE_001",
        }
    }

    # When first writing with add_to_maf=True
    write_variants_to_variant_file(variant_row, tn_pairs, True, variant_file)

    # Then
    # Check that the file now exists and contains the header + the variant row
    assert variant_file.exists(), "Variant file should be created."
    contents = variant_file.read_text().strip().split("\n")
    assert len(contents) == 2, "Should have a header and one variant line."
    assert (
        contents[0]
        == "TUMOUR\tNORMAL\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tAlternate_Allele\tAction"
    )
    expected_line = (
        "TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001\tTP53\t17\t7579472\t7579472\tC\tT\tADD"
    )
    assert contents[1] == expected_line, "Variant line should match expected."

    # When writing the same variant again
    write_variants_to_variant_file(variant_row, tn_pairs, True, variant_file)

    # Then
    # There should still be only one variant line (no duplicates)
    contents = variant_file.read_text().strip().split("\n")
    assert len(contents) == 2, "No duplicate variants should be appended."
    assert contents[1] == expected_line, "Duplicate line should not have been added."

    # When writing the same variant with add_to_maf=False (remove)
    write_variants_to_variant_file(variant_row, tn_pairs, False, variant_file)

    # Then
    # A remove action for the same variant should now be appended
    # since it differs by the 'Action' field.
    contents = variant_file.read_text().strip().split("\n")
    assert len(contents) == 3, "Should now have a remove action line appended."
    expected_remove_line = (
        "TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001\tTP53\t17\t7579472\t7579472\tC\tT\tREMOVE"
    )
    assert (
        contents[2] == expected_remove_line
    ), "The remove action line should match expected."


# Tests for process_true_positives


def test_process_true_positives_germline_classification(tmp_path, caplog):
    caplog.set_level(logging.INFO)

    # Given:
    # Create a mock mpileup file with:
    # - One TRUE_POSITIVE variant (in TUMOUR_SAMPLE_001)
    # - The same variant present as FALSE_NEGATIVE in multiple other pairs
    # - Normals with ALT reads above threshold
    mpileup_data = """# Comment line
Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumour_Seq_Allele2\tTumor_Sample_Barcode\tAlt_Count\tTot_Count\tStatus
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_001\t10\t100\tTRUE_POSITIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_001\t0\t50\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_002\t5\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_002\t4\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_003\t5\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_003\t4\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_004\t5\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_004\t4\t100\tFALSE_NEGATIVE
"""
    mpileup_file = tmp_path / "mpileup_test.tsv"
    mpileup_file.write_text(mpileup_data)

    # Create a tumour-normal pairs file that includes the samples:
    # We have T1-N1, T2-N2, T3-N3, T4-N4 pairs.
    # The TRUE_POSITIVE variant is in T1. The FALSE_NEGATIVE ones are in T2, T3, T4.
    # Normal samples have ALT reads >3 to qualify as germline.
    tn_pairs_data = """TUMOUR\tNORMAL
TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001
TUMOUR_SAMPLE_002\tNORMAL_SAMPLE_002
TUMOUR_SAMPLE_003\tNORMAL_SAMPLE_003
TUMOUR_SAMPLE_004\tNORMAL_SAMPLE_004
"""
    tn_pairs_file = tmp_path / "tn_pairs_test.tsv"
    tn_pairs_file.write_text(tn_pairs_data)

    # The output variant file
    variant_file = tmp_path / "variants_output.tsv"

    # Parameters:
    # Setting a low min_germline_tn_pairs to 3 means that if 3 or more pairs are flagged as FALSE_NEGATIVE,
    # we start considering germline classification. We have T2, T3, T4 as false negatives (3 pairs).
    # min_alt_norm_reads=3 means normals with >3 ALT reads are considered germline.
    # In this data, NORMAL_SAMPLE_002, NORMAL_SAMPLE_003, and NORMAL_SAMPLE_004 have 4 ALT reads each.
    min_germline_tn_pairs = 3
    min_alt_norm_reads = 3

    # When:
    process_true_positives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_germline_tn_pairs=min_germline_tn_pairs,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
    )

    # Then:
    # We expect that the TRUE_POSITIVE variant for TP53 has been reclassified as germline and thus should be marked as REMOVE.
    assert variant_file.exists()
    contents = variant_file.read_text().strip().split("\n")

    # Check the header and the REMOVE line
    assert len(contents) == 2, "There should be a header plus one variant line."
    assert contents[0].startswith(
        "TUMOUR\tNORMAL\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tAlternate_Allele\tAction"
    )
    expected_remove_line = (
        "TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001\tTP53\t17\t7579472\t7579472\tC\tT\tREMOVE"
    )
    assert (
        contents[1] == expected_remove_line
    ), f"Expected remove line:\n{expected_remove_line}\nbut got:\n{contents[1]}"

    # Check logs for expected messages
    assert "Variant flagged as germline" in caplog.text


def test_process_true_positives_not_germline_due_to_insufficient_pairs(
    tmp_path, caplog
):
    caplog.set_level(logging.INFO)

    # Given:
    # TRUE_POSITIVE in TUMOUR_SAMPLE_001
    # FALSE_NEGATIVE in TUMOUR_SAMPLE_002 and NORMAL_SAMPLE_002
    mpileup_data = """# Comment
Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumour_Seq_Allele2\tTumor_Sample_Barcode\tAlt_Count\tTot_Count\tStatus
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_001\t10\t100\tTRUE_POSITIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_001\t0\t50\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_002\t5\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_002\t2\t100\tFALSE_NEGATIVE
"""
    mpileup_file = tmp_path / "mpileup_insufficient_pairs.tsv"
    mpileup_file.write_text(mpileup_data)

    tn_pairs_data = """TUMOUR\tNORMAL
TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001
TUMOUR_SAMPLE_002\tNORMAL_SAMPLE_002
"""
    tn_pairs_file = tmp_path / "tn_pairs_insufficient_pairs.tsv"
    tn_pairs_file.write_text(tn_pairs_data)

    variant_file = tmp_path / "variants_output_insufficient_pairs.tsv"

    # Set thresholds so germline classification is not possible
    min_germline_tn_pairs = 3
    min_alt_norm_reads = 3

    # When
    process_true_positives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_germline_tn_pairs=min_germline_tn_pairs,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
    )

    # Then
    assert not variant_file.exists(), "No germline variants should be flagged."
    assert "Variant likely somatic based on initial criteria." in caplog.text


def test_process_true_positives_not_germline_due_to_norm_reads(tmp_path, caplog):
    caplog.set_level(logging.INFO)

    # Given:
    # TRUE_POSITIVE in TUMOUR_SAMPLE_001
    # FALSE_NEGATIVE in T2, T3, T4 and their normals, meeting min_germline_tn_pairs
    # But normals do not have >3 ALT reads, so no germline classification.
    mpileup_data = """# Comment
Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumour_Seq_Allele2\tTumor_Sample_Barcode\tAlt_Count\tTot_Count\tStatus
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_001\t10\t100\tTRUE_POSITIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_001\t0\t50\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_002\t5\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_002\t2\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_003\t5\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_003\t2\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_004\t5\t100\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_004\t2\t100\tFALSE_NEGATIVE
"""
    mpileup_file = tmp_path / "mpileup_norm_reads_fail.tsv"
    mpileup_file.write_text(mpileup_data)

    tn_pairs_data = """TUMOUR\tNORMAL
TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001
TUMOUR_SAMPLE_002\tNORMAL_SAMPLE_002
TUMOUR_SAMPLE_003\tNORMAL_SAMPLE_003
TUMOUR_SAMPLE_004\tNORMAL_SAMPLE_004
"""
    tn_pairs_file = tmp_path / "tn_pairs_norm_reads_fail.tsv"
    tn_pairs_file.write_text(tn_pairs_data)

    variant_file = tmp_path / "variants_output_norm_reads_fail.tsv"

    # Parameters
    min_germline_tn_pairs = 3
    min_alt_norm_reads = 3

    # When
    process_true_positives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_germline_tn_pairs=min_germline_tn_pairs,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
    )

    # Then
    assert not variant_file.exists(), "No germline variant should be output."


def test_process_true_positives_expected_mutations_in_output(
    mpileup_file, tn_pairs_file, tmp_path, caplog
):
    caplog.set_level(logging.DEBUG)

    # Define thresholds that we know from prior knowledge will classify at least one known variant as germline.
    # From the mock file, we know BRAF meets the criteria under these thresholds:
    min_germline_tn_pairs = 3
    min_alt_norm_reads = 3

    variant_file = tmp_path / "variants_output.tsv"

    # Run the process_true_positives function on the fixture data
    process_true_positives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_germline_tn_pairs=min_germline_tn_pairs,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
    )

    # After processing, we expect the known germline variant to appear in the variant file.
    expected_variant_line = (
        "TUMOUR_SAMPLE_D\tNORMAL_SAMPLE_D\tBRAF\t7\t140453136\t140453136\tT\tA\tREMOVE"
    )

    # Read the variant file and check if the expected line is present
    assert variant_file.exists(), "Variant file should be created."
    contents = variant_file.read_text().strip().split("\n")
    # Header line is always first
    assert (
        contents[0]
        == "TUMOUR\tNORMAL\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tAlternate_Allele\tAction"
    )
    assert (
        expected_variant_line in contents
    ), f"Expected germline variant line not found: {expected_variant_line}"

    # Optionally check logs to ensure correct classification messages appeared
    assert (
        "Variant flagged as germline" in caplog.text
    ), "Log should indicate the variant was flagged as germline."


# Tests for process_false_negatives
def test_process_false_negatives_adds_variant_when_criteria_met(
    mpileup_file, tn_pairs_file, tmp_path, caplog
):
    caplog.set_level(logging.INFO)

    # Given thresholds that allow a known false negative variant to be added
    # For example, assume that a variant in TUMOUR_SAMPLE_001 is FALSE_NEGATIVE with Alt_Count=10 in tumour
    # and its matched normal has Alt_Count=0, and we set min_alt_tum_reads=5, min_alt_norm_reads=3
    min_alt_tum_reads = 5
    min_alt_norm_reads = 3

    variant_file = tmp_path / "variants_output.tsv"

    # Create an empty set of germline variants to simulate no pre-existing germline exclusions
    germline_variant_ids = set()

    # When
    process_false_negatives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_alt_tum_reads=min_alt_tum_reads,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
        germline_variant_ids=germline_variant_ids,
    )

    # Then
    assert variant_file.exists(), "Variant file should be created."
    contents = variant_file.read_text().strip().split("\n")
    # Check the header
    assert (
        contents[0]
        == "TUMOUR\tNORMAL\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tAlternate_Allele\tAction"
    )
    # Check that at least one expected variant line is present
    # This line should match a known variant from your fixture data that meets the criteria
    # Example (adjust as per your actual data):
    expected_line = (
        "TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001\tTP53\tchr17\t7579472\t7579472\tG\tA\tADD"
    )
    assert (
        expected_line in contents
    ), f"Expected variant line not found: {expected_line}"

    # Check logs for confirmation
    assert (
        "Identified false negative variant suitable for addition to MAF" in caplog.text
    )


def test_process_false_negatives_no_add_when_criteria_not_met(
    mpileup_file, tn_pairs_file, tmp_path, caplog
):
    caplog.set_level(logging.INFO)

    # Set thresholds so that no variants meet the criteria
    # For example, set very high tumour ALT reads threshold and very low normal ALT reads threshold
    min_alt_tum_reads = 100000  # very high, no tumour should meet this
    min_alt_norm_reads = 0  # normal must have <0 is impossible, so no variants added
    variant_file = tmp_path / "variants_output.tsv"

    # Create an empty set of germline variants to simulate no pre-existing germline exclusions
    germline_variant_ids = set()

    # When
    process_false_negatives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_alt_tum_reads=min_alt_tum_reads,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
        germline_variant_ids=germline_variant_ids,
    )

    # Then
    if variant_file.exists():
        contents = variant_file.read_text().strip().split("\n")
        # If the file was created at all, it should only contain a header
        assert (
            len(contents) == 1
        ), "No variants should be added when criteria are not met."
    else:
        # File may not even be created if nothing is added
        pass

    # Check logs
    assert (
        "No false negative variants found. Nothing to add." in caplog.text
        or "No variants met the criteria" in caplog.text
    ), "Logs should indicate no variants were added."


def test_process_false_negatives_only_tumour_samples_added(
    mpileup_file, tn_pairs_file, tmp_path, caplog
):
    caplog.set_level(logging.INFO)

    # Choose thresholds that might tempt adding normal samples if not filtered out
    # For instance, if a normal sample had ALT_Count that looked like it met criteria, but we must ensure
    # that only tumour samples are processed.
    min_alt_tum_reads = 1
    min_alt_norm_reads = 5
    variant_file = tmp_path / "variants_output.tsv"

    # Create an empty set of germline variants to simulate no pre-existing germline exclusions
    germline_variant_ids = set()

    # When
    process_false_negatives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_alt_tum_reads=min_alt_tum_reads,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
        germline_variant_ids=germline_variant_ids,
    )

    # Then
    # Verify that all added variants (if any) belong to known tumour samples.
    if variant_file.exists():
        contents = variant_file.read_text().strip().split("\n")
        # First line is header
        for line in contents[1:]:
            fields = line.split("\t")
            tumour_sample = fields[0]
            # Check that this sample is listed as a TUMOUR in tn_pairs_file
            tn_pairs_df = tsv_to_df(tn_pairs_file)
            assert (
                tumour_sample in tn_pairs_df["TUMOUR"].values
            ), f"Sample {tumour_sample} is not a tumour sample."
    else:
        # If no variants added, test still passes since no normal samples were added
        pass


def test_prevent_germline_readdition(tmp_path, caplog):
    """
    Test that variants flagged as germline are not re-added as false-negative variants using mock data.
    """
    caplog.set_level(logging.DEBUG)

    # Step 1: Create mock mpileup and tumour-normal pairs data
    mpileup_data = """Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumour_Seq_Allele2\tTumor_Sample_Barcode\tAlt_Count\tTot_Count\tStatus
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_001\t10\t100\tTRUE_POSITIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_001\t0\t50\tTRUE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_002\t5\t50\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_002\t4\t50\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tTUMOUR_SAMPLE_003\t20\t50\tFALSE_NEGATIVE
TP53\t17\t7579472\t7579472\tC\tT\tNORMAL_SAMPLE_003\t1\t50\tFALSE_NEGATIVE
"""

    tn_pairs_data = """TUMOUR\tNORMAL
TUMOUR_SAMPLE_001\tNORMAL_SAMPLE_001
TUMOUR_SAMPLE_002\tNORMAL_SAMPLE_002
TUMOUR_SAMPLE_003\tNORMAL_SAMPLE_003
"""

    # Write mock data to temporary files
    mpileup_file = tmp_path / "mock_mpileup.tsv"
    tn_pairs_file = tmp_path / "mock_tn_pairs.tsv"
    variant_file = tmp_path / "mock_variant_output.tsv"
    mpileup_file.write_text(mpileup_data)
    tn_pairs_file.write_text(tn_pairs_data)

    # Step 2: Process true positives to flag a variant as germline (REMOVE)
    min_germline_tn_pairs = 1  # Low threshold to ensure germline detection
    min_alt_norm_reads = 3  # Normals with ALT_Count >3 are flagged

    # Process true positives and get the germline variants set
    germline_variant_ids = process_true_positives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_germline_tn_pairs=min_germline_tn_pairs,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
    )

    logging.debug(
        f"Germline variant IDs after processing true positives: {germline_variant_ids}"
    )

    # Assert that the germline variant was flagged for removal
    assert variant_file.exists(), "Variant file should be created."
    with open(variant_file, "r") as vf:
        contents = vf.readlines()
    assert len(contents) > 1, "The variant file should have a REMOVE entry."
    assert "REMOVE" in contents[-1], "The last action should be REMOVE."

    # Step 3: Process false negatives and ensure no re-addition of germline variants
    min_alt_tum_reads = 1  # Low threshold to ensure addition criteria are met
    min_alt_norm_reads = 3  # Low threshold to ensure addition criteria are met

    # Process false negatives
    process_false_negatives(
        mpileup_file=mpileup_file,
        tn_pairs_file=tn_pairs_file,
        min_alt_tum_reads=min_alt_tum_reads,
        min_alt_norm_reads=min_alt_norm_reads,
        variant_file=variant_file,
        germline_variant_ids=germline_variant_ids,
    )

    # Re-read the variant file after false-negative processing
    with open(variant_file, "r") as vf:
        updated_contents = vf.readlines()

    # Ensure no ADD actions are created for germline variants
    assert all(
        "ADD" not in line for line in updated_contents if "TP53" in line
    ), "Variants flagged as germline should not be re-added to the variant file."

    # Log validation
    assert (
        "Skipping variant" in caplog.text
    ), "Log should indicate that a germline variant was skipped during false-negative processing."
