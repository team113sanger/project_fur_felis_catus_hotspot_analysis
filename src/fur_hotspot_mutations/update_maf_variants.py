import argparse
import logging
import sys
from pathlib import Path
import typing as t

import pandas as pd

from utils.maf_utils import is_maf_format

__version__ = "1.0.0"


def setup_logging() -> None:
    """
    Set up logging for the script with DEBUG level by default.
    """
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logging.debug("Logging is set up with DEBUG level.")


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments
    """
    logging.debug("Setting up argument parser.")
    parser = argparse.ArgumentParser(
        description="This script takes a MAF file as input, adds/removes variants "
        "according to a variant file, and then writes the updated variants "
        "to a new MAF file."
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    parser.add_argument(
        "-i", "--input_maf", type=Path, required=True, help="Path to the input MAF file"
    )

    parser.add_argument(
        "-v",
        "--variant_file",
        type=Path,
        required=True,
        help="Path to the variant file containing variants to add or remove from the input MAF",
    )

    parser.add_argument(
        "-m",
        "--mpileup_file",
        type=Path,
        required=True,
        help="Path to the mpileup file contaning mpileup statistics for the variants in the input MAF and variant file",
    )

    parser.add_argument(
        "-o",
        "--output_maf",
        type=Path,
        required=True,
        help="Path to the output MAF file containing the updated variants",
    )

    args = parser.parse_args()
    logging.debug("Arguments parsed successfully:")
    logging.debug(f"  input_maf: {args.input_maf}")
    logging.debug(f"  variant_file: {args.variant_file}")
    logging.debug(f"  mpileup_file: {args.mpileup_file}")
    logging.debug(f"  output_maf: {args.output_maf}")
    return args


# General utility functions


def _validate_file_exists(file_path: Path) -> None:
    """
    Validates that the specified file exists.

    Args:
        file_path (Path): The file path to validate.

    Raises:
        ValueError: If the file does not exist.
    """
    logging.debug(f"Validating existence of file: {file_path}")
    if not file_path.exists():
        raise ValueError(
            f"Error: {file_path} does not exist. Please check the file path."
        )


def _check_required_columns(
    df: pd.DataFrame, required_columns: t.List[str], file_label: str
) -> None:
    """
    Check that the given DataFrame has all of the expected columns.

    Args:
        df (pd.DataFrame): The DataFrame to check.
        required_columns (List[str]): A list of columns that must be present.
        file_label (str): A label indicating which file is being checked (for logging).

    Raises:
        ValueError: If the DataFrame is missing any required columns.
    """
    logging.info(f"Checking {file_label} DataFrame for required columns ...")
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"{file_label} file is missing required column: {col}")
    logging.info(f"Successfully checked {file_label} file for required columns.")


def _convert_chromosome_df_column_to_string(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the Chromosome column in a DataFrame to a string.

    Args:
        df (pd.DataFrame): The DataFrame whose Chromosome column should be converted.

    Returns:
        pd.DataFrame: The updated DataFrame with Chromosome as string.

    Raises:
        KeyError: If 'Chromosome' column is not present in the given dataframe.
    """
    logging.debug("Converting 'Chromosome' column to string type.")
    if "Chromosome" in df.columns:
        df["Chromosome"] = df["Chromosome"].astype(str)
        logging.debug("Successfully converted 'Chromosome' column to string.")
    else:
        raise KeyError("'Chromosome' column is not present in the given dataframe.")
    return df


# Functions for validating and reading the input MAF file
def validate_and_load_maf_file(maf_file: Path) -> pd.DataFrame:
    """
    Checks if a file has valid MAF format, then loads it into a DataFrame.

    Args:
        maf_file (Path): The path to the MAF file.

    Returns:
        pd.DataFrame: The loaded MAF DataFrame.

    Raises:
        ValueError: If the MAF file is not valid or does not exist.
    """
    _validate_file_exists(maf_file)
    logging.info(f"Loading MAF file from {maf_file} into a DataFrame ...")

    if is_maf_format(maf_file):
        maf_df = pd.read_csv(maf_file, sep="\t", comment="#")
        _convert_chromosome_df_column_to_string(maf_df)
        logging.debug(f"MAF DataFrame preview:\n{maf_df.head()}")
        logging.info(f"Successfully loaded {maf_file} into a DataFrame.")
        logging.debug(f"MAF DataFrame shape: {maf_df.shape}")
        return maf_df
    else:
        raise ValueError("Please provide a valid MAF file")


# Functions for validating and reading the variant file


def _check_variant_df_action_column(variant_df: pd.DataFrame) -> None:
    """
    Checks that all values present in the Action column of the variant DataFrame are valid.

    Args:
        variant_df (pd.DataFrame): The variant DataFrame.

    Raises:
        ValueError: If the variant file contains unrecognized actions.
    """
    logging.info(
        "Checking that the Action column of the variant DataFrame is valid ..."
    )

    expected_actions = {"ADD", "REMOVE"}
    unique_actions = variant_df["Action"].unique()
    logging.debug(f"Unique actions found in variant file: {unique_actions}")

    for action in unique_actions:
        if action not in expected_actions:
            raise ValueError(
                f"Variant file contains unrecognised action {action}. "
                "Expected only ADD or REMOVE. Please check input data"
            )

    logging.info("Successfully validated the Action column.")


def _validate_variant_df(variant_df: pd.DataFrame) -> None:
    """
    Performs a series of helper functions that validate the structure of the variant DataFrame.

    Args:
        variant_df (pd.DataFrame): The DataFrame representing variants to add or remove.

    Raises:
        ValueError: If the DataFrame fails any validation checks.
    """
    logging.info("Validating the variant DataFrame ...")

    # Required columns for variant file
    required_columns = [
        "TUMOUR",
        "NORMAL",
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Alternate_Allele",
        "Action",
    ]
    _check_required_columns(variant_df, required_columns, "variant")

    _convert_chromosome_df_column_to_string(variant_df)
    _check_variant_df_action_column(variant_df)

    logging.info("Successfully validated the variant DataFrame.")


def _extract_variant_actions_from_variant_df(
    variant_df: pd.DataFrame,
) -> t.Tuple[bool, bool]:
    """
    Checks for the presence of ADD and REMOVE variants in the variant DataFrame.

    Args:
        variant_df (pd.DataFrame): The variant DataFrame.

    Returns:
        Tuple[bool, bool]: (contains_add_variants, contains_remove_variants)
    """
    contains_add_variants = "ADD" in variant_df["Action"].values
    contains_remove_variants = "REMOVE" in variant_df["Action"].values

    logging.debug(f"Variant DataFrame contains ADD variants: {contains_add_variants}")
    logging.debug(
        f"Variant DataFrame contains REMOVE variants: {contains_remove_variants}"
    )

    return (contains_add_variants, contains_remove_variants)


def validate_and_load_variant_file(
    variant_file: Path,
) -> t.Tuple[pd.DataFrame, t.Tuple[bool, bool]]:
    """
    Validates and loads the variant file.

    Args:
        variant_file (Path): The path to the variant file.

    Returns:
        Tuple[pd.DataFrame, Tuple[bool, bool]]:
            - The loaded and validated variant DataFrame.
            - A tuple indicating whether ADD/REMOVE variants are present.
    """
    _validate_file_exists(variant_file)
    logging.info(f"Loading variant file from {variant_file} into a DataFrame ...")
    variant_df = pd.read_csv(variant_file, sep="\t", comment="#")

    _validate_variant_df(variant_df)

    variant_actions = _extract_variant_actions_from_variant_df(variant_df)
    logging.debug(f"variant_actions: {variant_actions}")
    logging.debug(f"Variant DataFrame shape: {variant_df.shape}")
    return (variant_df, variant_actions)


# Functions for validating and reading the mpileup file
def validate_and_load_mpileup_file(mpileup_file: Path) -> pd.DataFrame:
    """
    Validates and loads the mpileup file into a DataFrame.

    Args:
        mpileup_file (Path): The path to the mpileup file.

    Returns:
        pd.DataFrame: The loaded mpileup DataFrame.

    Raises:
        ValueError: If required columns are missing.
    """
    _validate_file_exists(mpileup_file)
    logging.info(f"Loading mpileup file from {mpileup_file} into a DataFrame ...")
    mpileup_df = pd.read_csv(mpileup_file, sep="\t", comment="#")

    required_columns = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumour_Seq_Allele2",
        "Tumor_Sample_Barcode",
        "Alt_Count",
        "Tot_Count",
        "Alt_Perc",
    ]
    _check_required_columns(mpileup_df, required_columns, "mpileup")
    _convert_chromosome_df_column_to_string(mpileup_df)

    logging.debug(f"Mpileup DataFrame shape: {mpileup_df.shape}")
    return mpileup_df


# Functions for adding variants to the MAF file


def _extract_reference_maf_row(
    variant_row: pd.Series, maf_df: pd.DataFrame
) -> t.Dict[str, t.Any]:
    """
    Extract a 'reference' row from the MAF to serve as a template for the new variant.
    We only match on gene/position/alleles, ignoring the original Tumor_Sample_Barcode.

    Args:
        variant_row (pd.Series): A row from the variant file.
        maf_df (pd.DataFrame): The MAF DataFrame.

    Returns:
        dict: A dictionary representing a reference row from the MAF.

    Raises:
        ValueError: If no suitable reference row is found in the MAF.
    """
    logging.debug("Attempting to find a reference row in the MAF DataFrame.")
    matching_df = maf_df[
        (maf_df["Hugo_Symbol"] == variant_row["Hugo_Symbol"])
        & (maf_df["Chromosome"] == variant_row["Chromosome"])
        & (maf_df["Start_Position"] == variant_row["Start_Position"])
        & (maf_df["End_Position"] == variant_row["End_Position"])
        & (maf_df["Reference_Allele"] == variant_row["Reference_Allele"])
        & (maf_df["Tumor_Seq_Allele2"] == variant_row["Alternate_Allele"])
    ]

    if matching_df.empty:
        message = (
            "No suitable reference row found in the MAF for "
            f"{variant_row['Hugo_Symbol']} on {variant_row['Chromosome']} "
            f"{variant_row['Start_Position']}:{variant_row['End_Position']} "
            f"{variant_row['Reference_Allele']}>{variant_row['Alternate_Allele']}. "
            "Please ensure the MAF has at least one sample with this variant."
        )
        logging.error(message)
        raise ValueError(message)

    reference_row = matching_df.iloc[0].to_dict()
    logging.debug(
        f"Reference row found for {variant_row['Hugo_Symbol']}: {reference_row}"
    )

    return reference_row


def _calculate_new_maf_fields(
    variant_row: pd.Series, mpileup_file: Path
) -> t.Dict[str, t.Any]:
    """
    Use information from the mpileup value to calculate new fields for the given variant.

    Args:
        variant_row (pd.Series): A row from the variant file.
        mpileup_file (Path): The path to the mpileup file.

    Returns:
        dict: A dictionary of sample-specific data fields (counts, depths, etc.)

    Raises:
        ValueError: If tumour or normal data is empty.
    """
    logging.debug(
        f"Calculating new MAF fields for variant row: {variant_row.to_dict()}"
    )
    mpileup_df = validate_and_load_mpileup_file(mpileup_file)

    # Build boolean masks for tumor and normal
    tumor_mask = (
        (mpileup_df["Tumor_Sample_Barcode"] == variant_row["TUMOUR"])
        & (mpileup_df["Hugo_Symbol"] == variant_row["Hugo_Symbol"])
        & (mpileup_df["Chromosome"] == variant_row["Chromosome"])
        & (mpileup_df["Start_Position"] == variant_row["Start_Position"])
        & (mpileup_df["End_Position"] == variant_row["End_Position"])
        & (mpileup_df["Reference_Allele"] == variant_row["Reference_Allele"])
        & (mpileup_df["Tumour_Seq_Allele2"] == variant_row["Alternate_Allele"])
    )

    normal_mask = (
        (mpileup_df["Tumor_Sample_Barcode"] == variant_row["NORMAL"])
        & (mpileup_df["Hugo_Symbol"] == variant_row["Hugo_Symbol"])
        & (mpileup_df["Chromosome"] == variant_row["Chromosome"])
        & (mpileup_df["Start_Position"] == variant_row["Start_Position"])
        & (mpileup_df["End_Position"] == variant_row["End_Position"])
        & (mpileup_df["Reference_Allele"] == variant_row["Reference_Allele"])
        & (mpileup_df["Tumour_Seq_Allele2"] == variant_row["Alternate_Allele"])
    )

    tumour_variant_df = mpileup_df[tumor_mask]
    normal_variant_df = mpileup_df[normal_mask]

    # Check that both the tumour and normal samples have data
    if tumour_variant_df.empty:
        msg = "Tumour data is empty. Cannot create new MAF fields. Please check input data."
        logging.error(msg)
        raise ValueError(msg)
    if normal_variant_df.empty:
        msg = "Normal data is empty. Cannot create new MAF fields. Please check input data."
        logging.error(msg)
        raise ValueError(msg)

    # Calculate new MAF fields
    t_alt_count = tumour_variant_df["Alt_Count"].iloc[0]
    t_depth = tumour_variant_df["Tot_Count"].iloc[0]
    t_ref_count = t_depth - t_alt_count
    vaf_tum = tumour_variant_df["Alt_Perc"].iloc[0]

    n_alt_count = normal_variant_df["Alt_Count"].iloc[0]
    n_depth = normal_variant_df["Tot_Count"].iloc[0]
    n_ref_count = n_depth - n_alt_count
    vaf_norm = normal_variant_df["Alt_Perc"].iloc[0]

    sample_data = {
        "VAF_tum": vaf_tum,
        "t_ref_count": t_ref_count,
        "t_alt_count": t_alt_count,
        "t_depth": t_depth,
        "VAF_norm": vaf_norm,
        "n_ref_count": n_ref_count,
        "n_alt_count": n_alt_count,
        "n_depth": n_depth,
    }
    logging.debug(f"Calculated sample data for new variant: {sample_data}")
    return sample_data


def _construct_new_maf_row(
    reference_maf_row: t.Dict[str, t.Any], variant_row: pd.Series, mpileup_file: Path
) -> t.Dict[str, t.Any]:
    """
    Construct a new MAF row for a given variant using a reference MAF row and mpileup data.

    Args:
        reference_maf_row (dict): A dictionary representing a reference MAF row.
        variant_row (pd.Series): A row from the variant file.
        mpileup_file (Path): Path to the mpileup file.

    Returns:
        dict: A dictionary representing the newly constructed MAF row.
    """
    logging.debug("Constructing a new MAF row from reference row and variant data.")
    sample_data = _calculate_new_maf_fields(variant_row, mpileup_file)

    new_row = reference_maf_row.copy()

    # Update the reference MAF row with sample-specific values
    new_row.update(
        {
            "Tumor_Sample_Barcode": variant_row["TUMOUR"],
            "TUMOUR": variant_row["TUMOUR"],
            "Reference_Allele": variant_row["Reference_Allele"],
            "Tumor_Seq_Allele2": variant_row["Alternate_Allele"],
            "VAF_tum": sample_data["VAF_tum"],
            "t_ref_count": sample_data["t_ref_count"],
            "t_alt_count": sample_data["t_alt_count"],
            "t_depth": sample_data["t_depth"],
            "NORMAL": variant_row["NORMAL"],
            "VAF_norm": sample_data["VAF_norm"],
            "n_ref_count": sample_data["n_ref_count"],
            "n_alt_count": sample_data["n_alt_count"],
            "n_depth": sample_data["n_depth"],
        }
    )
    logging.debug(f"New MAF row constructed: {new_row}")
    return new_row


def add_variants_to_maf_df(
    maf_df: pd.DataFrame, add_variant_df: pd.DataFrame, mpileup_file: Path
) -> pd.DataFrame:
    """
    Construct and add new MAF rows for each variant specified in `add_variant_df`.
    Mark newly added variants in a special column "Added_By_MAF_Updater".

    Args:
        maf_df (pd.DataFrame): Original MAF DataFrame.
        add_variant_df (pd.DataFrame): Variants to add.
        mpileup_file (Path): Path to the mpileup file.

    Returns:
        pd.DataFrame: Updated MAF DataFrame with newly added variants.
    """
    logging.debug(f"Preparing to add {len(add_variant_df)} new variants to the MAF.")
    # Ensure the original MAF has a column that indicates if the variant was added by this script
    if "Added_By_MAF_Updater" not in maf_df.columns:
        # Existing rows were NOT added by this script
        maf_df["Added_By_MAF_Updater"] = False

    new_maf_rows = []

    logging.info("Constructing new MAF rows for variants in variant file ...")
    for idx, row in add_variant_df.iterrows():
        variant_info = (
            f"{row['Hugo_Symbol']}:"
            f"{row['Start_Position']}-{row['End_Position']}:"
            f"{row['Reference_Allele']}>{row['Alternate_Allele']}"
        )
        logging.info(f"Processing variant #{idx} -> {variant_info} ...")

        reference_maf_row = _extract_reference_maf_row(row, maf_df)
        new_maf_row = _construct_new_maf_row(reference_maf_row, row, mpileup_file)
        new_maf_rows.append(new_maf_row)

    logging.info("Adding new rows to MAF file ...")
    new_maf_rows_df = pd.DataFrame(new_maf_rows)

    # Mark newly added rows
    new_maf_rows_df["Added_By_MAF_Updater"] = True

    # Combine the original MAF with newly added rows
    updated_maf_df = pd.concat([maf_df, new_maf_rows_df], ignore_index=True)
    logging.debug(
        f"Updated MAF DataFrame shape after adding variants: {updated_maf_df.shape}"
    )

    return updated_maf_df


# Functions for removing variants from the MAF file


def remove_variants_from_maf_df(
    maf_df: pd.DataFrame, remove_variant_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Removes variants specified in `remove_variant_df` from the given `maf_df`,
    taking into account the sample name.

    Args:
        maf_df (pd.DataFrame): Original MAF DataFrame.
        remove_variant_df (pd.DataFrame): Variants to remove.

    Returns:
        pd.DataFrame: Updated MAF DataFrame with specified variants removed.
    """
    logging.debug(
        f"Preparing to remove {len(remove_variant_df)} variants from the MAF."
    )
    logging.info("Removing specified variants from MAF DataFrame...")

    # Convert remove_variant_df to a set of tuples for faster membership checks
    remove_set = set(
        zip(
            remove_variant_df["Hugo_Symbol"],
            remove_variant_df["Chromosome"],
            remove_variant_df["Start_Position"],
            remove_variant_df["End_Position"],
            remove_variant_df["Reference_Allele"],
            remove_variant_df["Alternate_Allele"],
            remove_variant_df["TUMOUR"],
        )
    )
    logging.debug(f"Number of variants to remove: {len(remove_set)}")

    to_remove_mask = maf_df.apply(
        lambda row: (
            (
                row["Hugo_Symbol"],
                row["Chromosome"],
                row["Start_Position"],
                row["End_Position"],
                row["Reference_Allele"],
                row["Tumor_Seq_Allele2"],
                row["Tumor_Sample_Barcode"],
            )
        )
        in remove_set,
        axis=1,
    )

    num_removed = to_remove_mask.sum()
    logging.info(f"Removing {num_removed} variants from MAF DataFrame.")

    updated_maf_df = maf_df[~to_remove_mask].reset_index(drop=True)
    logging.debug(f"Updated MAF DataFrame shape after removal: {updated_maf_df.shape}")
    return updated_maf_df


# Functions for processing variants in the variant file


def add_and_remove_variants_from_maf(
    variant_df: pd.DataFrame, maf_df: pd.DataFrame, mpileup_file: Path
) -> pd.DataFrame:
    """
    Remove and then add variants to the MAF DataFrame.

    Steps:
        1) Split 'variant_df' into remove-variants and add-variants
        2) Remove variants from 'maf_df'
        3) Add variants to 'maf_df'
        4) Return the updated 'maf_df'

    Args:
        variant_df (pd.DataFrame): DataFrame of variants containing ADD/REMOVE actions.
        maf_df (pd.DataFrame): Original MAF DataFrame.
        mpileup_file (Path): Path to the mpileup file.

    Returns:
        pd.DataFrame: Updated MAF DataFrame after removal and addition of variants.
    """
    logging.info("Performing BOTH add and remove actions.")
    add_variant_df = variant_df[variant_df["Action"] == "ADD"]
    remove_variant_df = variant_df[variant_df["Action"] == "REMOVE"]

    logging.debug(
        f"Number of REMOVE variants: {len(remove_variant_df)}; Number of ADD variants: {len(add_variant_df)}"
    )

    updated_maf_df = remove_variants_from_maf_df(maf_df, remove_variant_df)
    updated_maf_df = add_variants_to_maf_df(
        updated_maf_df, add_variant_df, mpileup_file
    )

    return updated_maf_df


def only_add_variants_to_maf(
    variant_df: pd.DataFrame, maf_df: pd.DataFrame, mpileup_file: Path
) -> pd.DataFrame:
    """
    Only add variants specified in the 'variant_df' (Action == 'ADD') to the MAF DataFrame.

    Args:
        variant_df (pd.DataFrame): DataFrame of variants with 'ADD' action.
        maf_df (pd.DataFrame): Original MAF DataFrame.
        mpileup_file (Path): Path to the mpileup file.

    Returns:
        pd.DataFrame: Updated MAF DataFrame with newly added variants.
    """
    logging.info("Performing ONLY add actions.")
    add_variant_df = variant_df[variant_df["Action"] == "ADD"]
    logging.debug(f"Number of variants to ADD: {len(add_variant_df)}")
    updated_maf_df = add_variants_to_maf_df(maf_df, add_variant_df, mpileup_file)
    return updated_maf_df


def only_remove_variants_from_maf(
    variant_df: pd.DataFrame, maf_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Only remove variants specified in the 'variant_df' (Action == 'REMOVE') from the MAF DataFrame.

    Args:
        variant_df (pd.DataFrame): DataFrame of variants with 'REMOVE' action.
        maf_df (pd.DataFrame): Original MAF DataFrame.

    Returns:
        pd.DataFrame: Updated MAF DataFrame with removed variants.
    """
    logging.info("Performing ONLY remove actions.")
    remove_variant_df = variant_df[variant_df["Action"] == "REMOVE"]
    logging.debug(f"Number of variants to REMOVE: {len(remove_variant_df)}")
    updated_maf_df = remove_variants_from_maf_df(maf_df, remove_variant_df)
    return updated_maf_df


def handle_no_variant_actions(
    variant_df: pd.DataFrame, maf_df: pd.DataFrame, mpileup_file: Path
) -> pd.DataFrame:
    """
    Handle the scenario where no variants need to be added or removed.

    Args:
        variant_df (pd.DataFrame): DataFrame of variants.
        maf_df (pd.DataFrame): Original MAF DataFrame.
        mpileup_file (Path): Path to the mpileup file.

    Returns:
        pd.DataFrame: Unchanged MAF DataFrame.
    """
    logging.info("No variants to add or remove. Returning input MAF unchanged.")
    return maf_df


def _process_variants_based_on_actions(
    variant_df: pd.DataFrame,
    variant_actions: t.Tuple[bool, bool],
    maf_df: pd.DataFrame,
    mpileup_file: Path,
) -> pd.DataFrame:
    """
    Decide which workflow to run based on whether ADD or REMOVE variants are present.

    Args:
        variant_df (pd.DataFrame): DataFrame containing variants with ADD/REMOVE actions.
        variant_actions (Tuple[bool, bool]): A tuple (contains_add_variants, contains_remove_variants).
        maf_df (pd.DataFrame): The original MAF DataFrame.
        mpileup_file (Path): Path to the mpileup file.

    Returns:
        pd.DataFrame: The updated MAF DataFrame after performing the relevant actions.
    """
    logging.debug(f"Processing variant actions: {variant_actions}")
    processing_instructions = {
        (True, True): add_and_remove_variants_from_maf,
        (True, False): only_add_variants_to_maf,
        (False, True): only_remove_variants_from_maf,
        (False, False): handle_no_variant_actions,
    }

    processing_func = processing_instructions.get(variant_actions)
    updated_maf_df = processing_func(variant_df, maf_df, mpileup_file)
    return updated_maf_df


def process_variants(
    input_maf: Path, variant_file: Path, mpileup_file: Path, output_maf: Path
) -> None:
    """
    Main driver function to orchestrate:
      1) Reading & validating MAF
      2) Reading & validating variant file
      3) Determining variant actions (ADD/REMOVE)
      4) Updating the MAF accordingly
      5) Writing the updated MAF to disk

    Args:
        input_maf (Path): Path to the input MAF file.
        variant_file (Path): Path to the variant file.
        mpileup_file (Path): Path to the mpileup file.
        output_maf (Path): Path to write the updated MAF file.
    """
    logging.info("Starting processing of variants...")

    maf_df = validate_and_load_maf_file(input_maf)
    variant_df, variant_actions = validate_and_load_variant_file(variant_file)

    updated_maf_df = _process_variants_based_on_actions(
        variant_df, variant_actions, maf_df, mpileup_file
    )

    # Sort the resulting MAF by chromosomal position if needed
    updated_maf_df = updated_maf_df.sort_values(
        by=["Chromosome", "Start_Position", "End_Position"],
        ascending=[True, True, True],
    )
    logging.debug(
        "MAF DataFrame is now sorted by Chromosome, Start_Position, and End_Position."
    )

    logging.info(f"Writing updated MAF to {output_maf} ...")
    updated_maf_df.to_csv(output_maf, sep="\t", index=False)

    logging.info("Successfully finished processing variants.")


def main() -> int:
    """
    Main entry point for command-line execution.

    Returns:
        int: Return code (0 for success, non-zero for errors).
    """
    try:
        setup_logging()
        args = parse_arguments()
        process_variants(
            args.input_maf, args.variant_file, args.mpileup_file, args.output_maf
        )
        return 0
    except ValueError as ve:
        logging.error(f"ValueError encountered: {ve}")
        return 1
    except Exception as e:
        logging.exception(f"An unexpected error occurred: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
