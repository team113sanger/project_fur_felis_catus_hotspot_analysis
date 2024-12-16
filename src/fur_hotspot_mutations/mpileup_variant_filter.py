import argparse
import csv
import logging
from pathlib import Path
from typing import Union

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
        "-v",
        "--vaf",
        type=float,
        default=1,
        help="Threshold for the minimum variant allele frequency (VAF) for a given variant to be considered "
        "i.e. with the default value of 1, any variant showing < 1% VAF will not be used when making decisions on variant filtering.",
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
        "i.e. with the default value of 3, if >3 tumour-normal pairs are flagged as false negative then report the true positive variant as a false positive.",
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


def _load_and_validate_tn_pairs(tn_pairs_file: Path) -> pd.DataFrame:
    _validate_tn_pairs_file(tn_pairs_file)
    return tsv_to_df(tn_pairs_file)


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


def _validate_samples(mpileup_df: pd.DataFrame, tn_pairs_df: pd.DataFrame) -> dict:
    _check_tn_pairs(mpileup_df, tn_pairs_df)
    return _extract_tn_pairs_from_df(tn_pairs_df)


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


def _validate_mpileup_df_columns(df: pd.DataFrame):
    """
    Validates that the required columns are present in the DataFrame.

    Parameters:
        df (pd.DataFrame): DataFrame to validate.
        required_columns (list): List of required column names.

    Raises:
        ValueError: If any required column is missing.
    """
    required_columns = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Tumor_Sample_Barcode",
        "Status",
        "Alt_Count",
    ]

    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")


def _load_and_validate_mpileup(mpileup_file: Path) -> pd.DataFrame:
    _validate_mpileup_file(mpileup_file)
    mpileup_df = tsv_to_df(mpileup_file)
    _validate_mpileup_df_columns(mpileup_df)

    return mpileup_df


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


# Functions for writing variants to the output variant file
def _get_matched_normal(tumour_sample: str, tn_pairs: dict) -> str:
    normal_sample = None

    for pair_id, samples in tn_pairs.items():
        if samples["TUMOUR"] == tumour_sample:
            normal_sample = samples["NORMAL"]
            break

    if not normal_sample:
        raise ValueError(
            f"Error when extracting matched normal for {tumour_sample}. Either tumour sample not in T/N pairs or T/N pairs are empty."
        )

    return normal_sample


def _convert_chromosome_to_string(chromosome: Union[int, str]) -> str:
    # Chromosome cannot be empty, raise an error if so
    is_empty_string = isinstance(chromosome, str) and chromosome.strip() == ""
    if chromosome is None or is_empty_string:
        raise ValueError(
            "Value for chromosome is empty. Check input data. Requires a non-empty string or an integer."
        )

    return str(chromosome)


def _construct_variant_file_row(
    variant_row: pd.Series, tn_pairs: dict, add_to_maf: bool
) -> str:
    # Populate the values for the new variant file row using information from the mpileup row
    try:
        tumour = variant_row["Tumor_Sample_Barcode"]
        # Get the ID of the matched normal
        normal = _get_matched_normal(tumour, tn_pairs)
        hugo_symbol = variant_row["Hugo_Symbol"]
        # If needed, convert the chromosome value to a string
        chrom = _convert_chromosome_to_string(variant_row["Chromosome"])
        start_pos = variant_row["Start_Position"]
        end_pos = variant_row["End_Position"]
        ref_allele = variant_row["Reference_Allele"]
        alt_allele = variant_row["Tumour_Seq_Allele2"]
        action = "ADD" if add_to_maf else "REMOVE"
    except KeyError:
        raise ValueError(
            f"Variant row is missing one or more fields. Please check input data: {variant_row}"
        )

    # Construct the new variant file row
    new_variant_file_row = (
        "\t".join(
            [
                tumour,
                normal,
                hugo_symbol,
                chrom,
                str(start_pos),
                str(end_pos),
                ref_allele,
                alt_allele,
                action,
            ]
        )
        + "\n"
    )

    return new_variant_file_row


def _counstruct_variant_id(variant_row: pd.Series) -> str:
    """
    Construct a unique variant ID for a given variant.

    Args:
        variant_row (pd.Series): A row containing variant data

    Returns:
        str: The unique variant ID for this variant
    """

    variant_id = "|".join(
        [
            variant_row["Hugo_Symbol"],
            _convert_chromosome_to_string(variant_row["Chromosome"]),
            str(variant_row["Start_Position"]),
            str(variant_row["End_Position"]),
            variant_row["Reference_Allele"],
            variant_row["Tumour_Seq_Allele2"],
        ]
    )

    return variant_id


def write_variants_to_variant_file(
    variant_row: pd.Series, tn_pairs: dict, add_to_maf: bool, variant_file: Path
) -> None:
    # If the variant file does not exist or is empty, write the expected header to the variant file path
    if not variant_file.exists() or variant_file.stat().st_size == 0:
        with open(variant_file, "w") as file:
            header = "TUMOUR\tNORMAL\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tAlternate_Allele\tAction\n"
            file.write(header)

    new_row = _construct_variant_file_row(variant_row, tn_pairs, add_to_maf)

    with open(variant_file, "r") as file:
        existing_rows = set(file.readlines())
    if new_row not in existing_rows:
        with open(variant_file, "a") as file:
            file.write(new_row)


# Functions for processing true positive variants
def _filter_false_negative_pairs(
    variant_rows_df: pd.DataFrame, tn_pairs: dict
) -> pd.DataFrame:
    """
    Filters the DataFrame to include only tumor-normal pairs where both samples have a FALSE_NEGATIVE mutation.

    Parameters:
        variant_rows_df (pd.DataFrame): DataFrame containing rows for a specific variant across samples.
        tn_pairs (dict): Dictionary of tumor-normal pairs.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only tumor-normal pairs with FALSE_NEGATIVE status.
    """
    # logging.debug(f"Input variant_rows_df:\n{variant_rows_df}")
    filtered_rows = []

    for pair_id, samples in tn_pairs.items():
        tumour_sample_id = samples["TUMOUR"]
        normal_sample_id = samples["NORMAL"]

        # Retrieve rows for the tumor and normal samples
        tumour_row = variant_rows_df[
            variant_rows_df["Tumor_Sample_Barcode"] == tumour_sample_id
        ]
        normal_row = variant_rows_df[
            variant_rows_df["Tumor_Sample_Barcode"] == normal_sample_id
        ]

        # logging.debug(f"Tumour row: {tumour_row}, Normal row: {normal_row}")

        # Ensure both rows exist and both have a FALSE_NEGATIVE status
        if not tumour_row.empty and not normal_row.empty:
            if (
                tumour_row["Status"].values[0] == "FALSE_NEGATIVE"
                and normal_row["Status"].values[0] == "FALSE_NEGATIVE"
            ):
                # Add both rows to the filtered list
                filtered_rows.append(tumour_row.iloc[0])
                filtered_rows.append(normal_row.iloc[0])

    result_df = pd.DataFrame(filtered_rows)
    # logging.debug(f"Filtered false negative pairs:\n{result_df}")

    return result_df


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


def _is_germline_pair(
    pair_id: str,
    tumour_status: str,
    normal_status: str,
    tumour_vaf: float,
    normal_vaf: float,
    tumour_alt_count: int,
    normal_alt_count: int,
    min_vaf: float,
    min_alt_norm_reads: int,
) -> bool:
    """
    Determines whether the given tumour-normal pair meets *any* germline criterion,
    and logs the result. Returns True if germline criteria are satisfied, False otherwise.
    """

    # Criterion 1: Tumour TRUE_NEGATIVE, Normal FALSE_NEGATIVE with > min_alt_norm_reads
    if (
        tumour_status == "TRUE_NEGATIVE"
        and normal_status == "FALSE_NEGATIVE"
        and normal_vaf >= min_vaf
        and normal_alt_count >= min_alt_norm_reads
    ):
        logging.debug(
            f"{pair_id} flagged as germline with {tumour_alt_count} ALT reads >= {min_vaf}% VAF in tumour and {normal_alt_count} ALT reads >= {min_vaf}% VAF in normal."
        )
        return True

    # Criterion 2: Tumour FALSE_NEGATIVE, Normal FALSE_NEGATIVE with > min_alt_norm_reads
    if (
        tumour_status == "FALSE_NEGATIVE"
        and normal_status == "FALSE_NEGATIVE"
        and normal_vaf >= min_vaf
        and normal_alt_count >= min_alt_norm_reads
    ):
        logging.debug(
            f"{pair_id} flagged as germline with {tumour_alt_count} ALT reads >= {min_vaf}% VAF in tumour and {normal_alt_count} ALT reads >= {min_vaf}% VAF in normal."
        )
        return True

    # Criterion 3: Tumour TRUE_POSITIVE, Normal FALSE_NEGATIVE
    if (
        tumour_status == "TRUE_POSITIVE"
        and normal_status == "FALSE_NEGATIVE"
        and normal_vaf >= min_vaf
        and normal_alt_count >= min_alt_norm_reads
    ):
        logging.debug(
            f"{pair_id} flagged as germline with {tumour_alt_count} ALT reads >= {min_vaf}% VAF in tumour and {normal_alt_count} ALT reads >= {min_vaf}% VAF in normal."
        )
        return True

    # Variant does not satisfy germline criteria
    logging.debug(
        f"{pair_id} is likely somatic with {tumour_alt_count} ALT reads >= {min_vaf}% VAF in tumour and {normal_alt_count} ALT reads >= {min_vaf}% VAF in normal."
    )
    return False


def _count_germline_tn_pairs(
    variant_rows_df: pd.DataFrame,
    tn_pairs: dict,
    min_vaf: float,
    min_alt_norm_reads: int,
) -> int:
    """
    Counts the number of tumour-normal pairs that satisfy germline criteria (see docstrings).
    Logs which pairs were flagged or skipped.
    """
    logging.info("Counting germline tumour-normal pairs with this mutation ...")

    germline_pair_count = 0

    for pair_id, samples in tn_pairs.items():
        tumour_sample_id = samples["TUMOUR"]
        normal_sample_id = samples["NORMAL"]

        # Retrieve rows for the tumour and normal samples
        tumour_row = variant_rows_df[
            variant_rows_df["Tumor_Sample_Barcode"] == tumour_sample_id
        ]
        normal_row = variant_rows_df[
            variant_rows_df["Tumor_Sample_Barcode"] == normal_sample_id
        ]

        if not tumour_row.empty and not normal_row.empty:
            tumour_status = tumour_row["Status"].values[0]
            normal_status = normal_row["Status"].values[0]
            tumour_vaf = tumour_row["Alt_Perc"].values[0]
            normal_vaf = normal_row["Alt_Perc"].values[0]
            tumour_alt_count = tumour_row["Alt_Count"].values[0]
            normal_alt_count = normal_row["Alt_Count"].values[0]

            if _is_germline_pair(
                pair_id,
                tumour_status,
                normal_status,
                tumour_vaf,
                normal_vaf,
                tumour_alt_count,
                normal_alt_count,
                min_vaf,
                min_alt_norm_reads,
            ):
                germline_pair_count += 1

    logging.info(f"Counted {germline_pair_count} germline tumour-normal pairs.")
    return germline_pair_count


def _log_flagged_germline_pairs(
    matching_variant_rows_df: pd.DataFrame,
    tn_pairs: dict,
    min_vaf: float,
    min_alt_norm_reads: int,
) -> None:
    flagged_pairs = []

    for pair_id, samples in tn_pairs.items():
        normal_sample = samples["NORMAL"]
        normal_row = matching_variant_rows_df[
            matching_variant_rows_df["Tumor_Sample_Barcode"] == normal_sample
        ]

        if not normal_row.empty:
            normal_status = normal_row["Status"].values[0]
            normal_vaf = normal_row["Alt_Perc"].values[0]
            normal_alt_count = normal_row["Alt_Count"].values[0]

            # Check if this normal sample truly contributed to germline detection
            if (
                normal_status == "FALSE_NEGATIVE"
                and normal_vaf >= min_vaf
                and normal_alt_count >= min_alt_norm_reads
            ):
                flagged_pairs.append(pair_id)

    logging.debug(f"Flagged pairs: {flagged_pairs}")

    if flagged_pairs:
        logging.debug(
            f"Flagged tumour-normal pairs as germline: {', '.join(flagged_pairs)}"
        )
    else:
        logging.debug("No tumour-normal pairs flagged as germline for this variant.")


def _process_germline_variants(
    mpileup_df: pd.DataFrame,
    tn_pairs: dict,
    min_vaf: float,
    min_germline_tn_pairs: int,
    min_alt_norm_reads: int,
    variant_file: Path,
):
    true_pos_df = extract_true_positives_from_mpileup_df(mpileup_df)
    logging.info(
        "Iterating through true positive variants to check for potential germline (and removal from MAF) ..."
    )

    germline_rows = []
    germline_variant_ids = set()

    for _, row in true_pos_df.iterrows():
        logging.info(
            f"Processing {row['Hugo_Symbol']} variant at position {row['Chromosome']}:{row['Start_Position']} "
            f"in sample {row['Tumor_Sample_Barcode']} ..."
        )
        matching_variant_rows_df = _extract_matching_variant_rows(row, mpileup_df)
        germline_tn_pair_count = _count_germline_tn_pairs(
            matching_variant_rows_df, tn_pairs, min_vaf, min_alt_norm_reads
        )
        logging.debug(
            f"Detected {germline_tn_pair_count} tumour-normal pairs matching germline criteria."
        )

        if germline_tn_pair_count < min_germline_tn_pairs:
            logging.info("Variant likely somatic based on initial criteria.")
            continue

        logging.info(
            f"Identified potential germline variant ({germline_tn_pair_count} out of {len(tn_pairs)} "
            "tumour-normal pairs have a false negative variant). Performing further checks ..."
        )

        if germline_tn_pair_count >= min_germline_tn_pairs:
            _log_flagged_germline_pairs(
                matching_variant_rows_df, tn_pairs, min_vaf, min_alt_norm_reads
            )
            logging.info(
                f"Variant flagged as germline (normal samples in {germline_tn_pair_count} tumour-normal pairs "
                f"have >= {min_alt_norm_reads} ALT reads). Will be removed from MAF"
            )
            germline_rows.append(row)

            #
            variant_id = _counstruct_variant_id(row)
            germline_variant_ids.add(variant_id)

    for germline_row in germline_rows:
        write_variants_to_variant_file(
            variant_row=germline_row,
            tn_pairs=tn_pairs,
            add_to_maf=False,
            variant_file=variant_file,
        )

    return germline_variant_ids


def process_true_positives(
    mpileup_file: Path,
    tn_pairs_file: Path,
    min_vaf: float,
    min_germline_tn_pairs: int,
    min_alt_norm_reads: int,
    variant_file: Path,
):
    logging.info(
        f"Checking true positive variants in {str(mpileup_file)} for potential germline variants ..."
    )
    mpileup_df = _load_and_validate_mpileup(mpileup_file)
    tn_pairs_df = _load_and_validate_tn_pairs(tn_pairs_file)
    tn_pairs = _validate_samples(mpileup_df, tn_pairs_df)
    germline_variant_ids = _process_germline_variants(
        mpileup_df,
        tn_pairs,
        min_vaf,
        min_germline_tn_pairs,
        min_alt_norm_reads,
        variant_file,
    )

    return germline_variant_ids


# Functions for processing false negative variants
def _process_single_variant(
    variant: pd.Series,
    mpileup_df: pd.DataFrame,
    tn_pairs: dict,
    min_vaf: float,
    min_alt_tum_reads: int,
    min_alt_norm_reads: int,
    variant_file: Path,
    tumour_samples: set,
) -> bool:
    """
    Processes a single false negative variant to determine if it meets the criteria for addition to the MAF file.

    Parameters:
        variant (pd.Series): A row from the false negative variants DataFrame.
        mpileup_df (pd.DataFrame): The mpileup DataFrame containing variant information.
        tn_pairs (dict): Dictionary mapping tumour-normal pairs.
        min_alt_tum_reads (int): Minimum ALT reads threshold in tumour samples.
        min_alt_norm_reads (int): Maximum ALT reads threshold in normal samples.
        variant_file (Path): Path to the output variant file.
        tumour_samples (set): Set of tumour sample IDs for quick lookup.

    Returns:
        bool: True if the variant was added to the MAF file, False otherwise.

    Raises:
        ValueError: If variant data for the tumour-normal pair is missing.
    """
    logging.debug(
        f"Processing variant {variant['Hugo_Symbol']} at {variant['Chromosome']}:{variant['Start_Position']} in {variant['Tumor_Sample_Barcode']} ..."
    )

    tumour_sample = variant["Tumor_Sample_Barcode"]

    # Ensure we're only processing tumour samples
    if tumour_sample not in tumour_samples:
        logging.info(
            f"Skipping variant {variant['Hugo_Symbol']} at {variant['Chromosome']}:{variant['Start_Position']} in {variant['Tumor_Sample_Barcode']} "
            f"because it originates from a normal sample."
        )
        return False

    # Get the matched normal sample
    normal_sample = _get_matched_normal(tumour_sample, tn_pairs)

    # Define common filter conditions
    filter_conditions = (
        (mpileup_df["Hugo_Symbol"] == variant["Hugo_Symbol"])
        & (mpileup_df["Chromosome"] == variant["Chromosome"])
        & (mpileup_df["Start_Position"] == variant["Start_Position"])
        & (mpileup_df["End_Position"] == variant["End_Position"])
        & (mpileup_df["Reference_Allele"] == variant["Reference_Allele"])
        & (mpileup_df["Tumour_Seq_Allele2"] == variant["Tumour_Seq_Allele2"])
    )

    # Retrieve variant data for the tumour sample
    tum_variant = mpileup_df[
        (mpileup_df["Tumor_Sample_Barcode"] == tumour_sample) & filter_conditions
    ]

    # Retrieve variant data for the normal sample
    norm_variant = mpileup_df[
        (mpileup_df["Tumor_Sample_Barcode"] == normal_sample) & filter_conditions
    ]

    if tum_variant.empty or norm_variant.empty:
        raise ValueError(
            f"Missing variant data for tumour-normal pair: {tumour_sample}-{normal_sample}."
        )

    # Count the ALT reads in the tumour and the normal
    vaf_tumour = tum_variant["Alt_Perc"].values[0]
    vaf_normal = norm_variant["Alt_Perc"].values[0]
    alt_count_tumour = tum_variant["Alt_Count"].values[0]
    alt_count_normal = norm_variant["Alt_Count"].values[0]

    logging.debug(f"# ALT reads in tumour: {alt_count_tumour}")
    logging.debug(f"Tumour variant VAF: {vaf_tumour}")
    logging.debug(f"# of ALT reads in normal: {alt_count_normal}")
    logging.debug(f"Normal variant VAF: {vaf_normal}")

    # Apply the criteria
    if (
        vaf_tumour >= min_vaf
        # and vaf_normal >= min_vaf
        and alt_count_tumour >= min_alt_tum_reads
        and alt_count_normal <= min_alt_norm_reads
    ):
        logging.info(
            f"Identified false negative variant suitable for addition to MAF: "
            f"{variant['Hugo_Symbol']} at {variant['Chromosome']}:{variant['Start_Position']} "
            f"in tumour {tumour_sample} (# of tumour ALT reads >= {min_vaf}% VAF: {alt_count_tumour}, "
            f"# of normal ALT reads >= {min_vaf}% VAF: {alt_count_normal})."
        )
        write_variants_to_variant_file(
            variant_row=tum_variant.iloc[0],
            tn_pairs=tn_pairs,
            add_to_maf=True,
            variant_file=variant_file,
        )
        return True

    return False


def _should_skip_variant(variant, germline_variant_ids):
    variant_id = _counstruct_variant_id(variant)
    if variant_id in germline_variant_ids:
        logging.info(
            f"Skipping variant {variant['Hugo_Symbol']} at {variant['Chromosome']}:{variant['Start_Position']} in {variant['Tumor_Sample_Barcode']} "
            f"as it has been previously flagged as germline."
        )
        return True
    return False


def _process_variants(
    false_neg_df,
    mpileup_df,
    tn_pairs,
    min_vaf,
    min_alt_tum_reads,
    min_alt_norm_reads,
    variant_file,
    tumour_samples,
    germline_variant_ids,
):
    variants_added = 0
    for _, variant in false_neg_df.iterrows():
        if _should_skip_variant(variant, germline_variant_ids):
            continue
        if _process_single_variant(
            variant=variant,
            mpileup_df=mpileup_df,
            tn_pairs=tn_pairs,
            min_vaf=min_vaf,
            min_alt_tum_reads=min_alt_tum_reads,
            min_alt_norm_reads=min_alt_norm_reads,
            variant_file=variant_file,
            tumour_samples=tumour_samples,
        ):
            variants_added += 1
    return variants_added


def _process_fn_variants(
    mpileup_df: pd.DataFrame,
    tn_pairs: dict,
    min_vaf: float,
    min_alt_tum_reads: int,
    min_alt_norm_reads: int,
    variant_file: Path,
    germline_variant_ids: set,
) -> None:
    """
    Processes false negative variants to identify those suitable for addition to the MAF file.
    """
    logging.info("Processing false negative variants for potential addition to MAF...")

    # Extract all false negative variants
    false_neg_df = extract_false_negatives_from_mpileup_df(mpileup_df)

    if false_neg_df.empty:
        logging.info("No false negative variants found. Nothing to add.")
        return

    # Prepare necessary data
    tumour_samples = {pair["TUMOUR"] for pair in tn_pairs.values()}
    germline_variant_ids = germline_variant_ids or set()

    # Process each variant
    variants_added = _process_variants(
        false_neg_df,
        mpileup_df,
        tn_pairs,
        min_vaf,
        min_alt_tum_reads,
        min_alt_norm_reads,
        variant_file,
        tumour_samples,
        germline_variant_ids,
    )

    logging.info(
        f"{'No variants met the criteria for addition to MAF.' if variants_added == 0 else f'Total variants added to MAF: {variants_added}'}"
    )


def process_false_negatives(
    mpileup_file: Path,
    tn_pairs_file: Path,
    min_vaf: float,
    min_alt_tum_reads: int,
    min_alt_norm_reads: int,
    variant_file: Path,
    germline_variant_ids: set,
):
    """
    High-level function to load data and processes false negative tumour variants to check if they meet the criteria for being added to the MAF file
    Args:
        mpileup_file (Path): Path to the mpileup file
        tn_pairs_file (Path): Path to the tumour-normal pairs file
        min_alt_tum_reads (int): Minimum number of ALT reads in the tumour to flag a variant as likely somatic
        min_alt_norm_reads (int): Minimum number of ALT reads in the matched normal to flag a variant as likely germline
        variant_file (Path): Output file to which variants that meet the criteria will be written
    """
    logging.info(
        f"Checking false negative variants in {str(mpileup_file)} for possible addition to MAF file ..."
    )

    mpileup_df = _load_and_validate_mpileup(mpileup_file)
    tn_pairs_df = _load_and_validate_tn_pairs(tn_pairs_file)
    tn_pairs = _validate_samples(mpileup_df, tn_pairs_df)

    _process_fn_variants(
        mpileup_df,
        tn_pairs,
        min_vaf,
        min_alt_tum_reads,
        min_alt_norm_reads,
        variant_file,
        germline_variant_ids,
    )


def main():
    args = parse_arguments()
    mpileup_file = args.mpileup_file
    tn_pairs_file = args.tn_pairs_file
    variant_file = args.variant_file
    min_vaf = args.vaf
    min_alt_tum_reads = args.alt_tumour_reads
    min_alt_norm_reads = args.alt_normal_reads
    min_germline_tn_pairs = args.germline_tn_pairs
    logging_level = args.log_level

    # Setup logging
    setup_logging(logging_level)

    # Process true positive variants and check if they are likely germline
    germline_variants = process_true_positives(
        mpileup_file,
        tn_pairs_file,
        min_vaf,
        min_germline_tn_pairs,
        min_alt_norm_reads,
        variant_file,
    )

    process_false_negatives(
        mpileup_file,
        tn_pairs_file,
        min_vaf,
        min_alt_tum_reads,
        min_alt_norm_reads,
        variant_file,
        germline_variants,
    )


if __name__ == "__main__":
    main()
