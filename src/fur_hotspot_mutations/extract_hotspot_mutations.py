import argparse
import logging
from pathlib import Path

import pandas as pd

from utils.maf_utils import is_maf_format


def setup_logging():
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract all hotspot mutations (i.e. mutations occurring at a hotspot site in > 2 samples) from a given MAF file. Outputs hotspot mutations to a new MAF file."
    )
    parser.add_argument(
        "-i",
        metavar="INPUT_MAF",
        nargs=1,
        type=Path,
        required=True,
        help="Path to input MAF file to extract hotspot mutations from.",
    )
    parser.add_argument(
        "-o",
        metavar="OUTPUT_MAF",
        nargs=1,
        type=Path,
        required=True,
        help="Path to output MAF with hotspot mutations extracted from -i.",
    )

    return parser.parse_args()


def extract_hotspot_mutations(maf_file, output_maf=None, save_to_file=True):
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

    # Filter for mutations occurring in at least two samples
    hotspot_mutations = grouped[grouped["Tumor_Sample_Barcode"] >= 2]

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


def main():
    # Parse command line arguments
    args = parse_arguments()

    input_maf = args.i[0]
    output_maf = args.o[0]

    # Check if input file is in MAF format
    if is_maf_format(input_maf):
        extract_hotspot_mutations(input_maf, output_maf)
    else:
        raise ValueError(
            f"{str(input_maf)} is not a valid MAF file. Please check input data."
        )


if __name__ == "__main__":
    setup_logging()
    main()
