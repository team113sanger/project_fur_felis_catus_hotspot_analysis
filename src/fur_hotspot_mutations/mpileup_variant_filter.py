import typing as t
import argparse
import csv
import logging
from pathlib import Path
from typing import Union
import warnings

import pandas as pd

from utils.logging_utils import setup_logging
from utils import constants


# CONSTANTS
# This constant is referenced by functions in other modules (e.g. cli.py)
# so do not factor it out.
COMMAND_NAME = constants.COMMAND_NAME__MPILEUP_VARIANT_FILTER

# FUNCTIONS


def get_argparser(
    subparser: t.Optional[argparse._SubParsersAction] = None,
) -> argparse.ArgumentParser:
    """
    Either returns a new ArgumentParser instance or a subparser for the
    mpileup_variant_filter command.

    A subparser is preferred, with unspecified behavior preserved for backwards
    compatibility.
    """
    description = constants.DESCRIPTION__MPILEUP_VARIANT_FILTER
    if subparser is None:
        parser = argparse.ArgumentParser(description=description)
    else:
        short_help = constants.SHORT_HELP__MPILEUP_VARIANT_FILTER
        parser = subparser.add_parser(
            COMMAND_NAME, description=description, help=short_help
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
        help="Path to the TSV file containing tumour normal pairs. "
        "This file should have the following headers: TUMOUR\tNORMAL",
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
        help="Threshold for the minimum variant allele frequency (VAF) for a given variant to be considered. "
        "i.e. with the default value of 1, any variant showing < 1 %% VAF will not be used for variant filtering decisions.",
    )

    parser.add_argument(
        "-t",
        "--alt_tumour_reads",
        type=int,
        default=5,
        help="Threshold for the minimum number of ALT allele tumour reads needed to report a given false negative variant "
        "as true. i.e. with the default value of 5, any FN variant with >5 ALT tumour reads (and <= -n normal reads) "
        "will be reported.",
    )

    parser.add_argument(
        "-n",
        "--alt_normal_reads",
        type=int,
        default=3,
        help="Threshold for the minimum number of ALT allele normal reads to classify a potential false negative as germline. "
        "i.e. with the default value of 3, any FN variant with >3 ALT normal reads is discarded as germline.",
    )

    parser.add_argument(
        "-g",
        "--germline_tn_pairs",
        type=int,
        default=3,
        help="Threshold for the minimum number of tumour-normal pairs flagged as false negatives to classify a true-positive variant as germline. "
        "i.e. with the default value of 3, if >3 T/N pairs are flagged FN, report the variant as germline (REMOVE).",
    )

    parser.add_argument(
        "--log-level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level (default: INFO)",
    )

    return parser


# -----------------------------------------------------------------------------
#  FILE-VALIDATION & DATAFRAME-LOADING
# -----------------------------------------------------------------------------
def _validate_file_exists(tsv_file: Path):
    """
    Validates that the specified file exists.

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
    """
    try:
        logging.debug(f"Loading TSV file {tsv_file} into DataFrame ...")
        return pd.read_csv(tsv_file, sep="\t", comment="#")
    except Exception as e:
        raise ValueError(f"Error reading {tsv_file}: {e}")


def tsv_to_df(tsv_file: Path) -> pd.DataFrame:
    """
    Reads a TSV file into a pandas DataFrame after removing comment lines.
    """
    _validate_file_exists(tsv_file)
    content = _get_non_comment_lines(tsv_file)
    _validate_delimiter(content, tsv_file)
    df = _load_tsv_to_dataframe(tsv_file)
    logging.debug(f"Successfully loaded {len(df)} rows from {tsv_file}")
    return df


# -----------------------------------------------------------------------------
#  TN-PAIRS FILE PROCESSING
# -----------------------------------------------------------------------------
def _validate_tn_pairs_file(tn_pairs_file: Path):
    """
    Validates that the tumour-normal pairs file exists and is not empty.
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
    logging.info(f"Loading and validating tumour-normal pairs file: {tn_pairs_file}")
    _validate_tn_pairs_file(tn_pairs_file)
    return tsv_to_df(tn_pairs_file)


def _check_tn_pairs(mpileup_df: pd.DataFrame, tn_pairs_df: pd.DataFrame):
    """
    Checks if the tumor and normal samples in the T/N pair DataFrame match the mpileup DataFrame samples.
    """
    logging.info(
        "Checking sample ID overlap between mpileup dataframe and tumour-normal pairs dataframe ..."
    )

    mpileup_samples = set(mpileup_df["Tumor_Sample_Barcode"])
    tn_pairs_samples = set(tn_pairs_df["TUMOUR"]).union(tn_pairs_df["NORMAL"])

    missing_in_tn_pairs = mpileup_samples - tn_pairs_samples
    if missing_in_tn_pairs:
        raise ValueError(
            f"Missing samples in tumour-normal pair file: {','.join(sorted(missing_in_tn_pairs))})"
        )

    missing_in_mpileup = tn_pairs_samples - mpileup_samples
    if missing_in_mpileup:
        raise ValueError(
            f"Missing samples in mpileup file: {','.join(sorted(missing_in_mpileup))}"
        )

    logging.info(
        "Successfully checked sample IDs. Both DataFrames share the same samples."
    )


def _extract_tn_pairs_from_df(tn_pairs_df: pd.DataFrame) -> dict:
    logging.debug("Extracting tumour-normal pairs from DataFrame ...")
    tn_pairs = {}
    for idx, row in tn_pairs_df.iterrows():
        pair_id = f"{row['TUMOUR']}_{row['NORMAL']}"
        tn_pairs[pair_id] = {"TUMOUR": row["TUMOUR"], "NORMAL": row["NORMAL"]}
    logging.debug(f"Extracted {len(tn_pairs)} tumour-normal pairs.")
    return tn_pairs


def _validate_samples(mpileup_df: pd.DataFrame, tn_pairs_df: pd.DataFrame) -> dict:
    _check_tn_pairs(mpileup_df, tn_pairs_df)
    return _extract_tn_pairs_from_df(tn_pairs_df)


# -----------------------------------------------------------------------------
#  MPILEUP FILE PROCESSING
# -----------------------------------------------------------------------------
def _validate_mpileup_file(mpileup_file: Path):
    try:
        if not mpileup_file.exists() or mpileup_file.stat().st_size == 0:
            raise ValueError(
                "Mpileup file is empty or does not exist. Please check input data."
            )
    except Exception as e:
        logging.error(f"Error validating mpileup file: {e}")
        raise


def _validate_mpileup_df_columns(df: pd.DataFrame):
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
    logging.info(f"Loading and validating mpileup file: {mpileup_file}")
    _validate_mpileup_file(mpileup_file)
    mpileup_df = tsv_to_df(mpileup_file)
    _validate_mpileup_df_columns(mpileup_df)
    return mpileup_df


def _extract_matching_variant_rows(row: pd.Series, mpileup_df: pd.DataFrame):
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
    logging.debug(
        f"Found {len(matching_rows_df)} matching rows for variant {gene} at {chrom}:{start_pos}-{end_pos}"
    )
    return matching_rows_df


def extract_false_negatives_from_mpileup_df(mpileup_df: pd.DataFrame) -> pd.DataFrame:
    false_neg_df = mpileup_df[mpileup_df["Status"] == "FALSE_NEGATIVE"]
    if false_neg_df.empty:
        logging.info("No false negatives detected. Skipping ...")
    else:
        logging.info(f"Detected {false_neg_df.shape[0]} false negative variants.")
    return false_neg_df


def extract_true_positives_from_mpileup_df(mpileup_df: pd.DataFrame) -> pd.DataFrame:
    true_pos_df = mpileup_df[mpileup_df["Status"] == "TRUE_POSITIVE"]
    if true_pos_df.empty:
        logging.info("No true positives detected. Skipping ...")
    else:
        logging.info(f"Detected {true_pos_df.shape[0]} true positive variants.")
    return true_pos_df


# -----------------------------------------------------------------------------
#  WRITING VARIANTS TO OUTPUT
# -----------------------------------------------------------------------------
def _get_matched_normal(tumour_sample: str, tn_pairs: dict) -> str:
    normal_sample = None
    for pair_id, samples in tn_pairs.items():
        if samples["TUMOUR"] == tumour_sample:
            normal_sample = samples["NORMAL"]
            break
    if not normal_sample:
        raise ValueError(
            f"Error when extracting matched normal for {tumour_sample}. "
            "Either tumour sample not in T/N pairs or T/N pairs are empty."
        )
    return normal_sample


def _convert_chromosome_to_string(chromosome: Union[int, str]) -> str:
    if chromosome is None or (isinstance(chromosome, str) and not chromosome.strip()):
        raise ValueError(
            "Value for chromosome is empty. Check input data. Requires a non-empty string or integer."
        )
    return str(chromosome)


def _construct_variant_file_row(
    variant_row: pd.Series, tn_pairs: dict, add_to_maf: bool
) -> str:
    try:
        tumour = variant_row["Tumor_Sample_Barcode"]
        normal = _get_matched_normal(tumour, tn_pairs)
        hugo_symbol = variant_row["Hugo_Symbol"]
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
    logging.debug(f"Constructed variant file row: {new_variant_file_row.strip()}")
    return new_variant_file_row


def _counstruct_variant_id(variant_row: pd.Series) -> str:
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
    if not variant_file.exists() or variant_file.stat().st_size == 0:
        with open(variant_file, "w") as file:
            header = "TUMOUR\tNORMAL\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tAlternate_Allele\tAction\n"
            file.write(header)
            logging.debug(f"Created variant file {variant_file} with header.")

    new_row = _construct_variant_file_row(variant_row, tn_pairs, add_to_maf)

    with open(variant_file, "r") as file:
        existing_rows = set(file.readlines())
    if new_row not in existing_rows:
        with open(variant_file, "a") as file:
            file.write(new_row)
        logging.info(
            f"{'Added' if add_to_maf else 'Removed'} variant row to {variant_file}: {new_row.strip()}"
        )


# -----------------------------------------------------------------------------
#  GERMLINE LOGIC (TRUE POSITIVES)
# -----------------------------------------------------------------------------
def _normal_meets_germline_criteria(
    normal_status: str,
    normal_vaf: float,
    normal_alt_count: int,
    min_vaf: float,
    min_alt_norm_reads: int,
) -> tuple[bool, list[str]]:
    """
    Checks if the normal sample meets the germline criteria:
    - normal_status == "FALSE_NEGATIVE"
    - normal_vaf >= min_vaf
    - normal_alt_count >= min_alt_norm_reads

    Returns (True, []) if all conditions are met, otherwise (False, [list_of_reasons]).
    """
    reasons = []

    if normal_status != "FALSE_NEGATIVE":
        reasons.append(f"Normal status {normal_status} != FALSE_NEGATIVE")
    if normal_vaf < min_vaf:
        reasons.append(f"Normal VAF {normal_vaf}% < {min_vaf}%")
    if normal_alt_count < min_alt_norm_reads:
        reasons.append(f"Normal ALT reads {normal_alt_count} < {min_alt_norm_reads}")

    meets_criteria = len(reasons) == 0
    return meets_criteria, reasons


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
    Returns True if the T/N pair meets any germline criterion, otherwise logs reasons and returns False.
    """

    # 1) Check if the normal sample meets germline criteria
    normal_ok, normal_reasons = _normal_meets_germline_criteria(
        normal_status,
        normal_vaf,
        normal_alt_count,
        min_vaf,
        min_alt_norm_reads,
    )

    # 2) Check tumour status
    tumour_ok_statuses = {"TRUE_NEGATIVE", "FALSE_NEGATIVE", "TRUE_POSITIVE"}
    tumour_ok = tumour_status in tumour_ok_statuses

    if normal_ok and tumour_ok:
        logging.debug(
            f"{pair_id} flagged as germline: "
            f"Tumour status={tumour_status}, Normal status={normal_status}, "
            f"Tumour VAF={tumour_vaf}%, Normal VAF={normal_vaf}% (>= {min_vaf}%), "
            f"Normal ALT reads={normal_alt_count} (>= {min_alt_norm_reads})."
        )
        return True

    # If not germline, gather reasons
    reasons = []
    if not normal_ok:
        reasons.extend(normal_reasons)
    if not tumour_ok:
        reasons.append(f"Tumour status {tumour_status} not in {tumour_ok_statuses}")

    logging.debug(
        f"{pair_id} is likely somatic: {', '.join(reasons) if reasons else 'No matching criteria'}"
    )
    return False


def _count_germline_tn_pairs(
    variant_rows_df: pd.DataFrame,
    tn_pairs: dict,
    min_vaf: float,
    min_alt_norm_reads: int,
) -> int:
    logging.info("Counting germline tumour-normal pairs with this mutation ...")
    germline_pair_count = 0

    for pair_id, samples in tn_pairs.items():
        tumour_sample_id = samples["TUMOUR"]
        normal_sample_id = samples["NORMAL"]

        tumour_row = variant_rows_df[
            variant_rows_df["Tumor_Sample_Barcode"] == tumour_sample_id
        ]
        normal_row = variant_rows_df[
            variant_rows_df["Tumor_Sample_Barcode"] == normal_sample_id
        ]

        if not tumour_row.empty and not normal_row.empty:
            tumour_status = tumour_row["Status"].values[0]
            normal_status = normal_row["Status"].values[0]
            tumour_vaf = tumour_row.get("Alt_Perc", 0).values[0]
            normal_vaf = normal_row.get("Alt_Perc", 0).values[0]
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
):
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

            if (
                normal_status == "FALSE_NEGATIVE"
                and normal_vaf >= min_vaf
                and normal_alt_count >= min_alt_norm_reads
            ):
                flagged_pairs.append(pair_id)

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

        if germline_tn_pair_count < min_germline_tn_pairs:
            shortfall = min_germline_tn_pairs - germline_tn_pair_count
            logging.info(
                f"Variant likely somatic. Germline pairs={germline_tn_pair_count} < min_germline_tn_pairs={min_germline_tn_pairs}. "
                f"Short by {shortfall} pair(s)."
            )
            continue

        logging.info(
            f"Identified potential germline variant ({germline_tn_pair_count} out of {len(tn_pairs)} "
            "tumour-normal pairs flagged as FN). Performing further checks ..."
        )

        _log_flagged_germline_pairs(
            matching_variant_rows_df, tn_pairs, min_vaf, min_alt_norm_reads
        )
        logging.info(
            f"Variant flagged as germline (normal samples in {germline_tn_pair_count} T/N pairs "
            f"have >= {min_alt_norm_reads} ALT reads). Will be removed from MAF"
        )
        germline_rows.append(row)
        germline_variant_ids.add(_counstruct_variant_id(row))

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
        f"Checking true positive variants in {mpileup_file} for potential germline variants ..."
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
    logging.info(
        f"Finished processing true positives. Identified {len(germline_variant_ids)} germline variant(s)."
    )
    return germline_variant_ids


# -----------------------------------------------------------------------------
#  FALSE NEGATIVES
# -----------------------------------------------------------------------------
def _meets_false_negative_criteria(
    vaf_tumour: float,
    alt_count_tumour: int,
    alt_count_normal: int,
    min_vaf: float,
    min_alt_tum_reads: int,
    min_alt_norm_reads: int,
) -> tuple[bool, list[str]]:
    """
    Checks whether a tumour variant meets the criteria for false negative reclassification.
    Returns (True, []) if it meets all criteria, otherwise (False, [reasons]).
    """
    reasons = []
    if vaf_tumour < min_vaf:
        reasons.append(f"tumour VAF {vaf_tumour}% < {min_vaf}%")
    if alt_count_tumour < min_alt_tum_reads:
        reasons.append(f"tumour ALT reads {alt_count_tumour} < {min_alt_tum_reads}")
    if alt_count_normal > min_alt_norm_reads:
        reasons.append(f"normal ALT reads {alt_count_normal} > {min_alt_norm_reads}")

    # If no reasons were accumulated, the criteria are met
    is_valid = len(reasons) == 0
    return is_valid, reasons


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
    logging.debug(
        f"Processing variant {variant['Hugo_Symbol']} at {variant['Chromosome']}:{variant['Start_Position']} "
        f"in {variant['Tumor_Sample_Barcode']} ..."
    )

    tumour_sample = variant["Tumor_Sample_Barcode"]
    if tumour_sample not in tumour_samples:
        logging.info(
            f"Skipping variant (from normal sample): {variant['Hugo_Symbol']} "
            f"at {variant['Chromosome']}:{variant['Start_Position']} "
            f"in {tumour_sample}."
        )
        return False

    normal_sample = _get_matched_normal(tumour_sample, tn_pairs)

    # Common filter conditions
    filter_conditions = (
        (mpileup_df["Hugo_Symbol"] == variant["Hugo_Symbol"])
        & (mpileup_df["Chromosome"] == variant["Chromosome"])
        & (mpileup_df["Start_Position"] == variant["Start_Position"])
        & (mpileup_df["End_Position"] == variant["End_Position"])
        & (mpileup_df["Reference_Allele"] == variant["Reference_Allele"])
        & (mpileup_df["Tumour_Seq_Allele2"] == variant["Tumour_Seq_Allele2"])
    )

    tum_variant = mpileup_df[
        (mpileup_df["Tumor_Sample_Barcode"] == tumour_sample) & filter_conditions
    ]
    norm_variant = mpileup_df[
        (mpileup_df["Tumor_Sample_Barcode"] == normal_sample) & filter_conditions
    ]

    if tum_variant.empty or norm_variant.empty:
        raise ValueError(
            f"Missing variant data for tumour-normal pair: {tumour_sample}-{normal_sample}."
        )

    vaf_tumour = tum_variant["Alt_Perc"].iloc[0]
    alt_count_tumour = tum_variant["Alt_Count"].iloc[0]
    alt_count_normal = norm_variant["Alt_Count"].iloc[0]

    meets_criteria, reasons = _meets_false_negative_criteria(
        vaf_tumour,
        alt_count_tumour,
        alt_count_normal,
        min_vaf,
        min_alt_tum_reads,
        min_alt_norm_reads,
    )

    if meets_criteria:
        logging.debug(
            f"Variant meets criteria for addition to MAF: Tumour VAF {vaf_tumour}% >= {min_vaf}%, "
            f"Tumour ALT reads {alt_count_tumour} >= {min_alt_tum_reads}, "
            f"Normal ALT reads {alt_count_normal} <= {min_alt_norm_reads}"
        )
        logging.info(
            f"Adding FN variant to MAF: {variant['Hugo_Symbol']} at {variant['Chromosome']}:{variant['Start_Position']} "
            f"in tumour {tumour_sample}."
        )
        write_variants_to_variant_file(
            variant_row=tum_variant.iloc[0],
            tn_pairs=tn_pairs,
            add_to_maf=True,
            variant_file=variant_file,
        )
        return True
    else:
        logging.info(
            f"Variant {variant['Hugo_Symbol']} at {variant['Chromosome']}:{variant['Start_Position']} "
            f"in tumour {tumour_sample} NOT added (reasons: {', '.join(reasons)})."
        )
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
    false_neg_df = extract_false_negatives_from_mpileup_df(mpileup_df)
    if false_neg_df.empty:
        logging.info("No false negative variants found. Nothing to add.")
        return

    tumour_samples = {pair["TUMOUR"] for pair in tn_pairs.values()}
    germline_variant_ids = germline_variant_ids or set()

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
    High-level function to load data and process false negative tumour variants to check if they meet the criteria for being added to the MAF file.
    """
    logging.info(
        f"Checking false negative variants in {mpileup_file} for possible addition to MAF file ..."
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


# -----------------------------------------------------------------------------
# MAIN PIPELINE
# -----------------------------------------------------------------------------
def main(args: argparse.Namespace):
    """
    Main pipeline to process mpileup data and identify questionable variants for
    potential addition/removal from a MAF file.
    """

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
    logging.info("Starting mpileup variant filtering script...")

    logging.info(
        f"User parameters: mpileup_file={mpileup_file}, tn_pairs_file={tn_pairs_file}, variant_file={variant_file}, "
        f"min_vaf={min_vaf}, alt_tumour_reads={min_alt_tum_reads}, alt_normal_reads={min_alt_norm_reads}, "
        f"germline_tn_pairs={min_germline_tn_pairs}"
    )

    # Step 1: Process true positive variants for potential germline reclassification
    germline_variants = process_true_positives(
        mpileup_file,
        tn_pairs_file,
        min_vaf,
        min_germline_tn_pairs,
        min_alt_norm_reads,
        variant_file,
    )

    # Step 2: Process false negative variants for possible addition to MAF
    process_false_negatives(
        mpileup_file,
        tn_pairs_file,
        min_vaf,
        min_alt_tum_reads,
        min_alt_norm_reads,
        variant_file,
        germline_variants,
    )

    logging.info("Completed mpileup variant filtering pipeline.")


if __name__ == "__main__":
    # The if-name-equals-main block is no-longer the preferred way to run scripts, as
    # we have moved to using a unified entry point for the CLI.
    #
    # However for backwards compatibility we show the user a deprecation warning.
    warnings.warn(
        (
            f"The script `{__file__}` should be run using the CLI program "
            f"running `{constants.PROGRAM_NAME}`."
        ),
        FutureWarning,
    )
    logging.info("Parsing command line arguments ...")
    parser = get_argparser(subparser=None)
    args = parser.parse_args()
    main(args)
