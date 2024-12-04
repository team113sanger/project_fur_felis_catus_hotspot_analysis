import argparse
import csv
import logging
from pathlib import Path

import pandas as pd


def setup_logging(level):
    """
    Configures the logging settings for the script.
    """
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def parse_arguments():
    """
    Parses command-line arguments provided by the user.

    Returns:
        argparse.Namespace: An object containing parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Uses mpileup data to identify questionable variant calls. Uses thresholds determined by the user to decide whether a given mutation should be included/excluded in the final callset."
    )

    parser.add_argument(
        "-m",
        "--mpileup_file",
        type=Path,
        required=True,
        help="Path to mpileup TSV file.",
    )

    parser.add_argument(
        "-p",
        "--tn_pairs_file",
        type=Path,
        required=True,
        help="Path to the TSV file containing tumour normal pairs. This file should have the following headers: TUMOUR\tNORMAL",
    )

    parser.add_argument(
        "-o",
        "--variant_file",
        type=Path,
        required=True,
        help="Path to the output variant file with variants to add/remove from a MAF file.",
    )

    parser.add_argument(
        "-t",
        "--alt_tumour_reads",
        type=int,
        default=5,
        help="Threshold for the minimum number of ALT allele tumour reads needed to report a given false negative variant as true "
        "i.e. with the default value of 5, any false negative variant showing >5 ALT tumour reads (and < -n normal reads; see below) will be reported.",
    )

    parser.add_argument(
        "-n",
        "--alt_normal_reads",
        type=int,
        default=3,
        help="Threshold for the minimum number of ALT allele normal reads to needed to classify a potential false negative as germline "
        "i.e. with the default value of 3, any false negative variant showing >3 ALT normal reads will be discarded.",
    )

    parser.add_argument(
        "-g",
        "--germline_tn_pairs",
        type=int,
        default=3,
        help="Threshold for the minimum number of tumour-normal pairs flagged as false negatives by mpileup to classify the true positive variants as germline "
        "i.e. with the default value of 5, if >5 tumour-normal pairs are flagged as false negative then report the true positive variant as a false positive.",
    )

    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level (default: INFO)",
    )

    return parser.parse_args()


# General functions
def _validate_file_exists(tsv_file: Path):
    """
    Validates that the specified file exists.

    Parameters:
        tsv_file (Path): The path to the TSV file.

    Raises:
        ValueError: If the file does not exist.
    """
    if not tsv_file.exists():
        raise ValueError(
            f"Error: {tsv_file} does not exist. Please check the file path."
        )


def _get_non_comment_lines(tsv_file: Path) -> str:
    """
    Reads the file and filters out lines starting with '#'.

    Parameters:
        tsv_file (Path): The path to the TSV file.

    Returns:
        str: The concatenated content of non-comment lines.

    Raises:
        ValueError: If the file is empty or contains only comment lines.
    """
    with open(tsv_file, "r") as file:
        lines = [line for line in file if not line.startswith("#")]
    if not lines:
        raise ValueError(
            f"Error: {tsv_file} is empty or contains only comments. Please check input data."
        )
    return "".join(lines)


def _validate_delimiter(content: str, tsv_file: Path):
    """
    Validates that the content of the file is tab-delimited.

    Parameters:
        content (str): The content of the file as a single string.
        tsv_file (Path): The path to the TSV file.

    Raises:
        ValueError: If the file is not tab-delimited or the delimiter cannot be determined.
    """
    try:
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(content)
        if dialect.delimiter != "\t":
            raise ValueError(
                f"Error: {tsv_file} is not tab-delimited. Please provide a tab-delimited file."
            )
    except csv.Error:
        raise ValueError(
            f"Error: Unable to determine the delimiter of {tsv_file}. Please check input data."
        )


def _load_tsv_to_dataframe(tsv_file: Path) -> pd.DataFrame:
    """
    Loads a tab-delimited TSV file into a pandas DataFrame.

    Parameters:
        tsv_file (Path): The path to the TSV file.

    Returns:
        pd.DataFrame: The DataFrame containing the TSV data.

    Raises:
        ValueError: If an error occurs while reading the file into a DataFrame.
    """
    try:
        return pd.read_csv(tsv_file, sep="\t", comment="#")
    except Exception as e:
        raise ValueError(f"Error reading {tsv_file}: {e}")


def tsv_to_df(tsv_file: Path) -> pd.DataFrame:
    """
    Reads a TSV file into a pandas DataFrame after removing comment lines.

    Parameters:
        tsv_file (Path): The path to the TSV file.

    Returns:
        pd.DataFrame: The DataFrame containing the TSV data.

    Raises:
        ValueError: If the file does not exist, is empty, or is not tab-delimited.
    """
    _validate_file_exists(tsv_file)
    content = _get_non_comment_lines(tsv_file)
    _validate_delimiter(content, tsv_file)
    return _load_tsv_to_dataframe(tsv_file)


# Functions for processing TN pairs file
def _validate_tn_pairs_file(tn_pairs_file: Path):
    """
    Validates that the tumour-normal pairs file exists and is not empty.

    Parameters:
        tn_pairs_file (Path): The path to the tumour-normal pairs file.

    Raises:
        ValueError: If the file does not exist or is empty.
    """
    try:
        if not tn_pairs_file.exists() or tn_pairs_file.stat().st_size == 0:
            raise ValueError(
                "Tumour-normal pairs file is empty or does not exist. Please check input data."
            )
    except Exception as e:
        logging.error(f"Error while validating tumour-normal pairs file: {e}")
        raise


def _check_tn_pairs(mpileup_df: pd.DataFrame, tn_pairs_df: pd.DataFrame):
    """
    Checks if the tumor and normal samples in the tumor-normal pair DataFrame match
    the samples found in the mpileup DataFrame and vice versa. Raises an error for any missing samples.

    Parameters:
        mpileup_df (pd.DataFrame): The mpileup DataFrame.
        tn_pairs_df (pd.DataFrame): The tumour-normal pairs DataFrame.

    Raises:
        ValueError: If any tumor or normal samples from the mpileup file are missing
                    in the tumor-normal pair file or vice versa.
    """
    logging.info(
        "Checking sample ID overlap between mpileup dataframe and tumour-normal pairs dataframe ..."
    )

    # Extract unique tumor and normal sample barcodes from the mpileup file
    mpileup_samples = set(mpileup_df["Tumor_Sample_Barcode"])

    # Load the tumor-normal pair file
    tn_pairs_samples = set(tn_pairs_df["TUMOUR"].values).union(
        set(tn_pairs_df["NORMAL"].values)
    )

    # Check for missing samples
    missing_samples_in_tn_pairs = mpileup_samples - tn_pairs_samples
    if missing_samples_in_tn_pairs:
        raise ValueError(
            f"Missing samples in tumour-normal pair file: {','.join(sorted(missing_samples_in_tn_pairs))})"
        )

    missing_samples_in_mpileup = tn_pairs_samples - mpileup_samples
    if missing_samples_in_mpileup:
        raise ValueError(
            f"Missing samples in mpileup file: {','.join(missing_samples_in_mpileup)}"
        )

    logging.info(
        "Successfully checked sample IDs. Both mpileup dataframe and tumour-normal pairs dataframe share the same samples."
    )


def _extract_tn_pairs_from_df(tn_pairs_df: pd.DataFrame) -> dict:
    """
    Extracts tumour-normal pairs from the given file and returns them as a dictionary.

    Parameters:
        tn_pairs_file (pd.DataFrame): The tumour-normal pairs DataFrame.

    Returns:
        dict: A dictionary containing the tumour-normal pairs.
              Structure: {pair_id: {'TUMOUR': tumour_sample_id, 'NORMAL': normal_sample_id}}
    """
    tn_pairs = {}
    for idx, row in tn_pairs_df.iterrows():
        pair_id = f"{row['TUMOUR']}_{row['NORMAL']}"
        tn_pairs[pair_id] = {"TUMOUR": row["TUMOUR"], "NORMAL": row["NORMAL"]}

    return tn_pairs


# Functions for processing mpileup file
def _validate_mpileup_file(mpileup_file: Path):
    """
    Validates that the mpileup file exists and is not empty.

    Parameters:
        mpileup_file (Path): The path to the mpileup file.

    Raises:
        ValueError: If the file does not exist or is empty.
    """
    try:
        if not mpileup_file.exists() or mpileup_file.stat().st_size == 0:
            raise ValueError(
                "Mpileup file is empty or does not exist. Please check input data."
            )
    except Exception as e:
        logging.error(f"Error validating mpileup file: {e}")
        raise


def _extract_matching_variant_rows(row: pd.Series, mpileup_df: pd.DataFrame):
    """
    Extracts rows from the mpileup DataFrame that match the given variant.

    Parameters:
        row (pd.Series): A series containing variant information.
        mpileup_df (pd.DataFrame): The mpileup DataFrame.

    Returns:
        pd.DataFrame: A DataFrame containing rows that match the variant.
    """
    gene = row["Hugo_Symbol"]
    chrom = row["Chromosome"]
    start_pos = row["Start_Position"]
    end_pos = row["End_Position"]
    ref_allele = row["Reference_Allele"]
    alt_allele = row["Tumour_Seq_Allele2"]

    matching_rows_df = mpileup_df[
        (mpileup_df["Hugo_Symbol"] == gene)
        & (mpileup_df["Chromosome"] == chrom)
        & (mpileup_df["Start_Position"] == start_pos)
        & (mpileup_df["End_Position"] == end_pos)
        & (mpileup_df["Reference_Allele"] == ref_allele)
        & (mpileup_df["Tumour_Seq_Allele2"] == alt_allele)
    ]

    return matching_rows_df


def extract_false_negatives_from_mpileup_df(mpileup_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts false negative variants from the mpileup DataFrame.

    Parameters:
        mpileup_df (pd.DataFrame): The mpileup DataFrame.

    Returns:
        pd.DataFrame: A DataFrame containing false negative variants.
    """
    false_neg_df = mpileup_df[mpileup_df["Status"] == "FALSE_NEGATIVE"]

    if false_neg_df.empty:
        logging.info("No false negatives detected. Skipping ...")
    else:
        logging.info(f"Detected {false_neg_df.shape[0]} false negative variants.")

    return false_neg_df


def extract_true_positives_from_mpileup_df(mpileup_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts true positive variants from the mpileup DataFrame.

    Parameters:
        mpileup_df (pd.DataFrame): The mpileup DataFrame.

    Returns:
        pd.DataFrame: A DataFrame containing true positive variants.
    """
    true_pos_df = mpileup_df[mpileup_df["Status"] == "TRUE_POSITIVE"]

    if true_pos_df.empty:
        logging.info("No true positives detected. Skipping ...")
    else:
        logging.info(f"Detected {true_pos_df.shape[0]} true positive variants.")

    return true_pos_df


# Functions for processing true positive variants
def _pair_is_present(
    sample_ids_in_df: set, tumour_sample_id: str, normal_sample_id: str
) -> bool:
    """
    Checks if both the tumour and normal sample IDs are present in the given set of sample IDs.

    Parameters:
        sample_ids_in_df (set): A set of sample IDs present in the DataFrame.
        tumour_sample_id (str): The tumour sample ID to check.
        normal_sample_id (str): The normal sample ID to check.

    Returns:
        bool: True if both sample IDs are present, False otherwise.
    """
    is_present = (
        tumour_sample_id in sample_ids_in_df and normal_sample_id in sample_ids_in_df
    )
    return is_present


def _count_tn_pairs(df: pd.DataFrame, tn_pairs_dict: dict) -> int:
    """
    Counts the number of tumour-normal pairs present in the given DataFrame.

    Parameters:
        df (pd.DataFrame): A DataFrame containing mpileup variant rows.
        tn_pairs_dict (dict): A dictionary containing tumour-normal pairs.
            Expected structure: {pair_id: {'TUMOUR': tumour_sample_id, 'NORMAL': normal_sample_id}}

    Returns:
        int: The number of tumour-normal pairs present in the DataFrame.
    """
    sample_ids_in_df = set(df["Tumor_Sample_Barcode"].unique())

    pair_count = 0

    for pair_id, samples in tn_pairs_dict.items():
        tumour_sample_id = samples["TUMOUR"]
        normal_sample_id = samples["NORMAL"]
        if _pair_is_present(sample_ids_in_df, tumour_sample_id, normal_sample_id):
            pair_count += 1

    return pair_count


def process_true_positives(
    mpileup_file: Path, tn_pairs_file: Path, min_germline_tn_pairs: int
):
    """
    Processes true positive variants from the mpileup file to identify potential germline variants.

    Parameters:
        mpileup_file (Path): The path to the mpileup file.
        tn_pairs_file (Path): The path to the tumour-normal pairs file.
        min_germline_tn_pairs (int): The minimum number of tumour-normal pairs flagged as false negatives
                                     to classify the true positive variant as germline.
    """
    # Load in mpileup file data
    logging.info(
        f"Checking true positive variants in {str(mpileup_file)} for potential germline variants ..."
    )
    _validate_mpileup_file(mpileup_file)
    mpileup_df = tsv_to_df(mpileup_file)

    # Load in tumour-normal pairs file
    _validate_tn_pairs_file(tn_pairs_file)
    tn_pairs_df = tsv_to_df(tn_pairs_file)

    # Ensure samples match between the mpileup file and the tn pairs file
    _check_tn_pairs(mpileup_df, tn_pairs_df)

    # Get a dictionary containing the tumour-normal pairs
    tn_pairs = _extract_tn_pairs_from_df(tn_pairs_df)

    # Extract true positive variants from mpileup file
    true_pos_df = extract_true_positives_from_mpileup_df(mpileup_df)

    # Check for potential germline variants in the true positives
    logging.info(
        "Iterating through true positive variants to check for potential germline ..."
    )
    for _, row in true_pos_df.iterrows():
        logging.info(
            f"Processing {row['Hugo_Symbol']} variant at position {row['Chromosome']}:{row['Start_Position']} in sample {row['Tumor_Sample_Barcode']} ..."
        )
        # Create a filtered dataframe containing only this variant across all samples
        matching_variant_rows_df = _extract_matching_variant_rows(row, mpileup_df)
        false_neg_variant_rows = matching_variant_rows_df[
            matching_variant_rows_df["Status"] == "FALSE_NEGATIVE"
        ]

        # Count the number of tumour-normal pairs with false negative mutations
        false_neg_tn_pair_count = _count_tn_pairs(false_neg_variant_rows, tn_pairs)
        logging.debug(
            f"Detected {false_neg_tn_pair_count} tumour-normal pairs with false negative mutations."
        )

        # Check if this variant satisfies our threshold
        if false_neg_tn_pair_count >= min_germline_tn_pairs:
            logging.info(
                f"Identified potential germline variant ({false_neg_tn_pair_count} out of {len(tn_pairs)} tumour-normal pairs have a false negative variant). Performing further checks ..."
            )

        else:
            logging.info("Successfully processed variant - likely not germline.")


def main():
    args = parse_arguments()
    mpileup_file = args.mpileup_file
    tn_pairs_file = args.tn_pairs_file
    # min_alt_tum_reads = args.alt_tumour_reads
    # min_alt_norm_reads = args.alt_normal_reads
    min_germline_tn_pairs = args.germline_tn_pairs
    logging_level = args.log_level

    # Setup logging
    setup_logging(logging_level)

    # Process true positive variants and check if they are likely germline
    process_true_positives(mpileup_file, tn_pairs_file, min_germline_tn_pairs)


if __name__ == "__main__":
    main()
