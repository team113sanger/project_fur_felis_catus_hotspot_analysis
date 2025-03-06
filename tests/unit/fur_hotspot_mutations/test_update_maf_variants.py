import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch
import tempfile

from fur_hotspot_mutations.update_maf_variants import (
    validate_and_load_maf_file,
    validate_and_load_variant_file,
    validate_and_load_mpileup_file,
    remove_variants_from_maf_df,
    add_variants_to_maf_df,
    _process_variants_based_on_actions,
    add_and_remove_variants_from_maf,
    only_add_variants_to_maf,
    only_remove_variants_from_maf,
    handle_no_variant_actions,
    process_variants,
)

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------


@pytest.fixture
def valid_maf_df():
    """
    MAF with three 'reference' rows:
      1) TP53 (TUM1)  -- so we can remove it
      2) EGFR (REF_EGFR)  -- reference row for eventually adding EGFR as TUM2
      3) KRAS (REF_KRAS)  -- reference row for eventually adding KRAS as TUM3
    """
    data = {
        "Hugo_Symbol": ["TP53", "EGFR", "KRAS"],
        "Chromosome": ["17", "7", "12"],
        "Start_Position": [7574000, 55249071, 25398284],
        "End_Position": [7574001, 55249072, 25398285],
        "Reference_Allele": ["C", "T", "G"],
        "Tumor_Seq_Allele2": ["G", "A", "T"],
        "Tumor_Sample_Barcode": ["TUM1", "REF_EGFR", "REF_KRAS"],
        "Center": ["DummyCenter", "DummyCenter", "DummyCenter"],
        "NCBI_Build": ["GRCh38", "GRCh38", "GRCh38"],
        "Variant_Classification": [
            "Missense_Mutation",
            "Missense_Mutation",
            "Missense_Mutation",
        ],
        "Variant_Type": ["SNP", "SNP", "SNP"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def valid_variants_df():
    """
    - One REMOVE row for TP53 (TUM1).
    - Two ADD rows: EGFR (TUM2) and KRAS (TUM3).
    """
    data = {
        "TUMOUR": ["TUM1", "TUM2", "TUM3"],
        "NORMAL": ["NORM1", "NORM2", "NORM3"],
        "Hugo_Symbol": ["TP53", "EGFR", "KRAS"],
        "Chromosome": ["17", "7", "12"],
        "Start_Position": [7574000, 55249071, 25398284],
        "End_Position": [7574001, 55249072, 25398285],
        "Reference_Allele": ["C", "T", "G"],
        "Alternate_Allele": ["G", "A", "T"],
        "Action": ["REMOVE", "ADD", "ADD"],
    }
    return pd.DataFrame(data)


@pytest.fixture
def valid_mpileup_df():
    """
    Mpileup lines for:
      - TP53 (TUM1 / NORM1)
      - EGFR (TUM2 / NORM2)
      - KRAS (TUM3 / NORM3)

    so that any "add" operation can find matching tumour/normal data
    for the new sample (TUM2, TUM3).
    We do NOT remove mpileup rows for 'REF_EGFR' or 'REF_KRAS'
    because we only match on gene/pos/alleles, not tumor sample.
    """
    data = {
        "Hugo_Symbol": [
            # TP53
            "TP53",
            "TP53",
            # EGFR
            "EGFR",
            "EGFR",
            # KRAS
            "KRAS",
            "KRAS",
        ],
        "Chromosome": [
            "17",
            "17",
            "7",
            "7",
            "12",
            "12",
        ],
        "Start_Position": [
            7574000,
            7574000,
            55249071,
            55249071,
            25398284,
            25398284,
        ],
        "End_Position": [
            7574001,
            7574001,
            55249072,
            55249072,
            25398285,
            25398285,
        ],
        "Reference_Allele": [
            "C",
            "C",
            "T",
            "T",
            "G",
            "G",
        ],
        "Tumour_Seq_Allele2": [
            "G",
            "G",
            "A",
            "A",
            "T",
            "T",
        ],
        "Tumor_Sample_Barcode": [
            "TUM1",
            "NORM1",
            "TUM2",
            "NORM2",
            "TUM3",
            "NORM3",
        ],
        "Alt_Count": [
            40,
            2,
            10,
            1,
            12,
            0,
        ],
        "Tot_Count": [
            100,
            100,
            100,
            100,
            200,
            200,
        ],
        "Alt_Perc": [
            40.0,
            2.0,
            10.0,
            1.0,
            6.0,
            0.0,
        ],
    }
    return pd.DataFrame(data)


# --------------------------------------------------------------------------
# Tests: Validation & Loading
# --------------------------------------------------------------------------


def test_validate_and_load_maf_file(tmp_path, valid_maf_df):
    """Test that validate_and_load_maf_file works with a valid MAF file."""
    # Write the fixture MAF DataFrame to disk
    maf_file = tmp_path / "test_valid_maf.maf"
    valid_maf_df.to_csv(maf_file, sep="\t", index=False)

    loaded_df = validate_and_load_maf_file(maf_file)

    # Assertions
    pd.testing.assert_frame_equal(loaded_df, valid_maf_df)


def test_validate_and_load_variant_file(tmp_path, valid_variants_df):
    """Test that validate_and_load_variant_file works with a valid variants file."""
    variant_file = tmp_path / "test_valid_variants.tsv"
    valid_variants_df.to_csv(variant_file, sep="\t", index=False)

    # This should not raise any errors
    loaded_df, actions = validate_and_load_variant_file(variant_file)

    pd.testing.assert_frame_equal(loaded_df, valid_variants_df)
    # Expected actions: (True, True) because we have "ADD" & "REMOVE"
    assert actions == (True, True)


def test_validate_and_load_mpileup_file(tmp_path, valid_mpileup_df):
    """Test that validate_and_load_mpileup_file works with a valid mpileup file."""
    mpileup_file = tmp_path / "test_valid_mpileup.tsv"
    valid_mpileup_df.to_csv(mpileup_file, sep="\t", index=False)

    # This should not raise any errors
    loaded_mpileup_df = validate_and_load_mpileup_file(mpileup_file)
    pd.testing.assert_frame_equal(loaded_mpileup_df, valid_mpileup_df)


# --------------------------------------------------------------------------
# Tests: Removing Variants
# --------------------------------------------------------------------------


def test_remove_variants_from_maf_df(valid_maf_df):
    """Test removing variants from MAF."""
    remove_data = {
        "TUMOUR": ["TUM1"],
        "Hugo_Symbol": ["TP53"],
        "Chromosome": ["17"],
        "Start_Position": [7574000],
        "End_Position": [7574001],
        "Reference_Allele": ["C"],
        "Alternate_Allele": ["G"],
    }
    remove_variant_df = pd.DataFrame(remove_data)

    updated_df = remove_variants_from_maf_df(valid_maf_df, remove_variant_df)

    # The row with TP53/TUM1 is removed => 2 rows left
    assert len(updated_df) == 2
    assert "TP53" not in updated_df["Hugo_Symbol"].values


# --------------------------------------------------------------------------
# Tests: Adding Variants
# --------------------------------------------------------------------------


def test_add_variants_to_maf_df(valid_maf_df, valid_mpileup_df, tmp_path):
    # MAF starts with 3 rows
    # We'll add a new row for a TP53 mutation in TUM1
    add_data = {
        "TUMOUR": ["TUM1"],
        "NORMAL": ["NORM1"],
        "Hugo_Symbol": ["TP53"],
        "Chromosome": ["17"],
        "Start_Position": [7574000],
        "End_Position": [7574001],
        "Reference_Allele": ["C"],
        "Alternate_Allele": ["G"],
        "Action": ["ADD"],
    }
    add_variant_df = pd.DataFrame(add_data)

    mpileup_file = tmp_path / "dummy_mpileup.tsv"
    valid_mpileup_df.to_csv(mpileup_file, sep="\t", index=False)

    # This adds a fourth row - the initial 3 + 1 new row
    new_df = add_variants_to_maf_df(valid_maf_df, add_variant_df, mpileup_file)
    assert len(new_df) == 4


# --------------------------------------------------------------------------
# Tests: Combined Add + Remove
# --------------------------------------------------------------------------


def test_add_and_remove_variants_from_maf(
    valid_maf_df, valid_variants_df, valid_mpileup_df, tmp_path
):
    """
    MAF starts with 3 rows: TP53 (TUM1), EGFR (REF_EGFR), KRAS (REF_KRAS).
    - Remove TP53 (TUM1) => 2 rows left (EGFR, KRAS reference).
    - Add EGFR (TUM2) and KRAS (TUM3) => 2 more rows, total = 4.
    """

    mpileup_file = tmp_path / "dummy_mpileup.tsv"
    valid_mpileup_df.to_csv(mpileup_file, sep="\t", index=False)

    updated_df = add_and_remove_variants_from_maf(
        valid_variants_df, valid_maf_df, mpileup_file
    )

    # Final row count = 4
    assert len(updated_df) == 4

    # The original row for TP53 (TUM1) is gone
    assert not (
        (updated_df["Hugo_Symbol"] == "TP53")
        & (updated_df["Tumor_Sample_Barcode"] == "TUM1")
    ).any(), "TP53 (TUM1) should have been removed."

    # The new EGFR (TUM2) and KRAS (TUM3) rows are present
    # (the reference rows for EGFR(REF_EGFR) and KRAS(REF_KRAS) still remain too)
    new_tumours = set(updated_df["Tumor_Sample_Barcode"])
    assert "TUM2" in new_tumours, "EGFR (TUM2) was not added."
    assert "TUM3" in new_tumours, "KRAS (TUM3) was not added."


# --------------------------------------------------------------------------
# Tests: Only Add
# --------------------------------------------------------------------------


def test_only_add_variants_to_maf(
    valid_maf_df, valid_variants_df, valid_mpileup_df, tmp_path
):
    """
    Test if only_add_variants_to_maf adds all 'ADD' variants (EGFR, KRAS)
    while ignoring the 'REMOVE' variant (TP53).

    Starting MAF has 3 rows:
      1) TP53 (TUM1)
      2) EGFR (REF_EGFR)
      3) KRAS (REF_KRAS)

    The 'valid_variants_df' has:
      - TP53 (TUM1) as REMOVE
      - EGFR (TUM2) and KRAS (TUM3) as ADD

    Because we are only adding variants, we expect:
      - The original 3 rows remain
      - 2 new rows (for EGFR(TUM2) and KRAS(TUM3)) are appended
      => total 5 rows in updated MAF.
    """

    # Write the mpileup fixture to a temporary file
    mpileup_file = tmp_path / "test_mpileup.tsv"
    valid_mpileup_df.to_csv(mpileup_file, sep="\t", index=False)

    # Call the function under test
    updated_df = only_add_variants_to_maf(
        variant_df=valid_variants_df, maf_df=valid_maf_df, mpileup_file=mpileup_file
    )

    # We started with 3 rows in the MAF, we add 2 variants => total 5
    assert len(updated_df) == len(valid_maf_df) + 2

    # Check that the REMOVE row (TP53, TUM1) is still present,
    # because 'only_add_variants_to_maf' does not remove anything.
    # The easiest check is to see if TUM1 is in the final Tumor_Sample_Barcode column:
    assert (
        "TUM1" in updated_df["Tumor_Sample_Barcode"].values
    ), "The original TP53 (TUM1) row should not be removed."

    # Check that new rows for TUM2 (EGFR) and TUM3 (KRAS) are appended
    # i.e., their Tumor_Sample_Barcodes are in the updated DataFrame
    new_tumours = set(updated_df["Tumor_Sample_Barcode"].values)
    assert "TUM2" in new_tumours, "EGFR (TUM2) should have been added."
    assert "TUM3" in new_tumours, "KRAS (TUM3) should have been added."


# --------------------------------------------------------------------------
# Tests: Only Remove
# --------------------------------------------------------------------------


def test_only_remove_variants_from_maf(
    valid_maf_df, valid_variants_df, valid_mpileup_df
):
    """
    Test if only_remove_variants_from_maf removes all 'REMOVE' variants and ignores 'ADD'.
    """
    # We have 1 'REMOVE' row in valid_variants_df. The function should only remove those.
    updated_df = only_remove_variants_from_maf(
        valid_variants_df, valid_maf_df, valid_mpileup_df
    )

    # The remove variant has TP53 with TUM1 => that row is removed
    # The row for "EGFR" in valid_maf_df is still there
    assert len(updated_df) == len(valid_maf_df) - 1


# --------------------------------------------------------------------------
# Tests: Handle No Variant Actions
# --------------------------------------------------------------------------


def test_handle_no_variant_actions(valid_maf_df):
    """
    If no variants have ADD or REMOVE, the returned MAF should be the same.
    """
    # Create a DataFrame with random columns but no "ADD"/"REMOVE"
    no_actions_df = pd.DataFrame(
        {
            "TUMOUR": ["TUM1"],
            "NORMAL": ["NORM1"],
            "Hugo_Symbol": ["TP53"],
            "Chromosome": ["17"],
            "Start_Position": [7574000],
            "End_Position": [7574001],
            "Reference_Allele": ["C"],
            "Alternate_Allele": ["G"],
            "Action": ["FOO"],  # No ADD or REMOVE
        }
    )

    updated_df = handle_no_variant_actions(
        no_actions_df, valid_maf_df, Path("dummy_mpileup.tsv")
    )

    pd.testing.assert_frame_equal(updated_df, valid_maf_df)


# --------------------------------------------------------------------------
# Tests: _process_variants_based_on_actions
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "actions, expected_func",
    [
        ((True, True), add_and_remove_variants_from_maf),
        ((True, False), only_add_variants_to_maf),
        ((False, True), only_remove_variants_from_maf),
        ((False, False), handle_no_variant_actions),
    ],
)
def test__process_variants_based_on_actions(
    valid_maf_df, valid_mpileup_df, actions, expected_func
):
    from fur_hotspot_mutations.update_maf_variants import (
        _process_variants_based_on_actions,
    )

    with patch(
        f"fur_hotspot_mutations.update_maf_variants.{expected_func.__name__}",
        return_value="dummy_return",
    ) as mocked_func:
        dummy_variants_df = pd.DataFrame()
        mpileup_file = Path("dummy_mpileup.tsv")

        retval = _process_variants_based_on_actions(
            dummy_variants_df, actions, valid_maf_df, mpileup_file
        )

        mocked_func.assert_called_once()
        assert retval == "dummy_return"


# --------------------------------------------------------------------------
# Tests: process_variants (End-to-end)
# --------------------------------------------------------------------------


def test_process_variants_end_to_end(
    tmp_path, valid_maf_df, valid_variants_df, valid_mpileup_df
):
    """
    An integration-like test for process_variants that:
      1) Writes sample input MAF, variants, mpileup to temp files
      2) Calls process_variants (removing TP53 TUM1, adding EGFR TUM2 & KRAS TUM3)
      3) Checks output MAF for expected changes
    """
    # Step 1) Write input files to disk
    maf_file = tmp_path / "input.maf"
    variants_file = tmp_path / "variants.tsv"
    mpileup_file = tmp_path / "mpileup.tsv"
    output_maf_file = tmp_path / "output.maf"

    valid_maf_df.to_csv(maf_file, sep="\t", index=False)
    valid_variants_df.to_csv(variants_file, sep="\t", index=False)
    valid_mpileup_df.to_csv(mpileup_file, sep="\t", index=False)

    # Step 2) Run the process_variants function
    process_variants(maf_file, variants_file, mpileup_file, output_maf_file)

    # Step 3) Load the output file and verify
    updated_df = pd.read_csv(output_maf_file, sep="\t")

    # Our MAF starts with 3 rows:
    #   1) TP53 (TUM1)
    #   2) EGFR (REF_EGFR)
    #   3) KRAS (REF_KRAS)
    #
    # By default (True, True) => add_and_remove_variants_from_maf
    #  * Removes TP53 (TUM1)
    #  * Adds EGFR (TUM2) and KRAS (TUM3)
    # So final row count => (start=3) - 1 + 2 = 4
    expected_count = len(valid_maf_df) - 1 + 2
    assert (
        len(updated_df) == expected_count
    ), f"Expected {expected_count} rows, got {len(updated_df)}."

    # Check that TP53 (TUM1) is removed
    assert not (
        (updated_df["Hugo_Symbol"] == "TP53")
        & (updated_df["Tumor_Sample_Barcode"] == "TUM1")
    ).any(), "TP53(TUM1) should have been removed but is still present."

    # Check that newly added EGFR (TUM2) and KRAS (TUM3) are present
    new_samples = set(updated_df["Tumor_Sample_Barcode"].values)
    assert "TUM2" in new_samples, "EGFR(TUM2) was not added to the MAF output."
    assert "TUM3" in new_samples, "KRAS(TUM3) was not added to the MAF output."

    # Check that original reference rows (EGFR(REF_EGFR) & KRAS(REF_KRAS)) remain
    assert (
        "REF_EGFR" in new_samples
    ), "Original EGFR(REF_EGFR) row should still be present."
    assert (
        "REF_KRAS" in new_samples
    ), "Original KRAS(REF_KRAS) row should still be present."
