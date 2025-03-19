import typing as t
import argparse
import logging
from pathlib import Path
import warnings

import pandas as pd

from utils.logging_utils import setup_logging
from utils.maf_utils import is_maf_format
from utils import constants

# CONSTANTS
# This constant is referenced by functions in other modules (e.g. cli.py)
# so do not factor it out.
COMMAND_NAME = constants.COMMAND_NAME__EXTRACT_HOTSPOT_MUTATIONS

# FUNCTIONS


def get_argparser(
    subparser: t.Optional[argparse._SubParsersAction] = None,
) -> argparse.ArgumentParser:
    """
    Either returns a new ArgumentParser instance or a subparser for the extract_hotspot_mutations command.

    A subparser is preferred, with unspecified behavior preserved for backwards
    compatibility.
    """
    description = constants.DESCRIPTION__EXTRACT_HOTSPOT_MUTATIONS
    if subparser is None:
        parser = argparse.ArgumentParser(description=description)
    else:
        short_help = constants.SHORT_HELP__EXTRACT_HOTSPOT_MUTATIONS
        parser = subparser.add_parser(
            COMMAND_NAME, description=description, help=short_help
        )

    parser.add_argument(
        "-i",
        metavar="INPUT_MAF",
        nargs=1,
        type=Path,
        required=True,
        help="Path to the input MAF file from which hotspot mutations will be extracted.",
    )
    parser.add_argument(
        "-m",
        metavar="MIN_SAMPLES",
        nargs=1,
        type=int,
        default=[2],
        help=(
            "Minimum number of samples required for a mutation to be classified as a hotspot mutation. "
            "Default value: 2."
        ),
    )
    parser.add_argument(
        "-o",
        metavar="OUTPUT_MAF",
        nargs=1,
        type=Path,
        required=True,
        help="Path to the output MAF file where extracted hotspot mutations will be saved.",
    )

    return parser


def extract_hotspot_mutations(
    maf_file, min_samples=2, output_maf=None, save_to_file=True
):
    """
    Extract mutations from a MAF file that meet the hotspot criteria:
    - Same gene
    - Same genomic position
    - Same mutation type
    - Occur in more than two samples

    Parameters:
    maf_file (str): Path to the input MAF file.
    output_maf (str, optional): Path to save the new MAF file containing hotspot mutations.
    save_to_file (bool): Whether to save the result to a file. Defaults to True.

    Returns:
    pd.DataFrame: DataFrame of the hotspot mutations.
    """
    logging.info(f"Extracting hotspot mutations from {str(maf_file)} ...")

    # Load the MAF file into a DataFrame
    maf_df = pd.read_csv(maf_file, sep="\t", comment="#", low_memory=False)

    # Ensure required columns are present
    required_columns = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode",
    ]
    for col in required_columns:
        if col not in maf_df.columns:
            raise ValueError(f"Required column '{col}' not found in the MAF file.")

    # Remove duplicate rows for the same mutation-sample pair
    maf_df = maf_df.drop_duplicates(
        subset=[
            "Hugo_Symbol",
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode",
        ]
    )

    # Group by mutation-specific attributes and count unique samples
    grouped = (
        maf_df.groupby(
            [
                "Hugo_Symbol",
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ]
        )
        .agg({"Tumor_Sample_Barcode": pd.Series.nunique})
        .reset_index()
    )

    # Filter for mutations occurring in at least the minimum number of samples
    logging.info(
        f"Identifying hotspot muations present in at least {min_samples} samples ..."
    )
    hotspot_mutations = grouped[grouped["Tumor_Sample_Barcode"] >= min_samples]

    # Merge with the original MAF to retain all details for hotspot mutations
    hotspot_maf = pd.merge(
        maf_df,
        hotspot_mutations[
            [
                "Hugo_Symbol",
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ]
        ],
        on=[
            "Hugo_Symbol",
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
        ],
        how="inner",
    )

    # Sort by chromosome and position
    hotspot_maf["Chromosome"] = hotspot_maf["Chromosome"].astype(str)
    hotspot_maf["Start_Position"] = hotspot_maf["Start_Position"].astype(int)
    hotspot_maf.sort_values(by=["Chromosome", "Start_Position"], inplace=True)

    # Save the hotspot mutations to a new MAF file (if required)
    if save_to_file and output_maf:
        hotspot_maf.to_csv(output_maf, sep="\t", index=False)
        logging.info(
            f"Output MAF file with hotspot mutations written to {str(output_maf)}"
        )

    return hotspot_maf


def main(args: argparse.Namespace) -> None:
    setup_logging()
    # Parse command line arguments
    input_maf = args.i[0]
    min_samples = args.m[0]
    output_maf = args.o[0]

    logging.debug(f"Input MAF: {input_maf}")
    logging.debug(f"Minimum number of samples: {min_samples}")
    logging.debug(f"Output MAF: {output_maf}")

    # Check if input file is in MAF format
    if is_maf_format(input_maf):
        extract_hotspot_mutations(
            input_maf, min_samples=min_samples, output_maf=output_maf, save_to_file=True
        )
    else:
        raise ValueError(
            f"{str(input_maf)} is not a valid MAF file. Please check input data."
        )


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
    parser = get_argparser(subparser=None)
    args = parser.parse_args()
    main(args)
