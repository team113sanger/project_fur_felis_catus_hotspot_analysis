#!/usr/bin/env python3

import sys
import os
import subprocess
import statistics
import re
import argparse
import logging

# --- Logging Setup ---


def setup_logging():
    """
    Configures the logging settings.
    Logs are output to stderr with a specific format and levels.
    """
    logger = logging.getLogger("VariantProcessor")
    logger.setLevel(logging.DEBUG)  # Set to DEBUG to capture all levels of logs

    # Create handlers
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)  # Capture all levels in the handler

    # Create formatter and add it to the handler
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)

    # Add handler to the logger
    logger.addHandler(handler)

    return logger


# --- Command-line Argument Setup


def parse_arguments():
    """
    Parses command-line arguments and returns them.
    """
    parser = argparse.ArgumentParser(
        description="Process MAF and BAM files to analyze genetic variants using samtools mpileup."
    )
    parser.add_argument(
        "bam_list_file", help="Path to the BAM list file containing paths to BAM files."
    )
    parser.add_argument(
        "maf_list_file", help="Path to the MAF list file containing paths to MAF files."
    )
    parser.add_argument("reference_fasta", help="Path to the reference FASTA file.")
    parser.add_argument(
        "map_quality",
        nargs="?",
        type=int,
        default=20,
        help="Mapping quality threshold (default: 20).",
    )
    parser.add_argument(
        "base_quality",
        nargs="?",
        type=int,
        default=20,
        help="Base quality threshold (default: 20).",
    )
    return parser.parse_args()


# --- File Validation and Reading ---


def check_required_files(file_paths, logger):
    """
    Checks if all required files exist. Exits the script with an error message if any file is missing.
    """
    logger.info("Checking for the existence of required files.")
    for file_path in file_paths:
        if not os.path.isfile(file_path):
            logger.error(f"File not found: {file_path}")
            sys.exit(f"Error: File {file_path} does not exist.\n")
        else:
            logger.debug(f"File exists: {file_path}")


def read_file_list(list_file, logger, file_type="MAF"):
    """
    Reads a list of file paths from the provided list file.
    """
    logger.info(f"Reading {file_type} list from file: {list_file}")
    try:
        with open(list_file, "r") as file:
            files = [line.strip() for line in file if line.strip()]
        logger.debug(f"{file_type} files found: {files}")
        return files
    except Exception as e:
        logger.error(f"Failed to read {file_type} list file {list_file}: {e}")
        sys.exit(f"Error: Failed to read {file_type} list file {list_file}: {e}\n")


# --- Variant Type Validation ---


def is_supported_variant_type(variant_type, logger):
    """
    Checks if the variant type is supported.
    Returns True if supported, else False.
    """
    supported_types = ["DEL", "INS", "DNP", "SNP", "TNP"]
    if variant_type not in supported_types:
        logger.debug(f"Unsupported variant type encountered: {variant_type}")
        return False
    return True


# --- Position Adjustment Functions ---


def adjust_ins_position(start_position, end_position):
    """
    Adjusts positions for insertion (INS) variants.
    """
    end_position = start_position
    return start_position, end_position


def adjust_del_position(start_position, end_position):
    """
    Adjusts positions for deletion (DEL) variants.
    """
    start_position -= 1
    end_position = start_position
    return start_position, end_position


def adjust_default_position(start_position, end_position):
    """
    Default position adjustment for variants that do not require specific adjustments.
    """
    return start_position, end_position


def adjust_positions(variant_type, start_position, end_position, logger):
    """
    Adjusts start and end positions based on variant type.
    Returns a tuple of (adjusted_start, adjusted_end) or None if variant type is unsupported.
    """
    ADJUSTMENT_FUNCTIONS = {"INS": adjust_ins_position, "DEL": adjust_del_position}
    try:
        adjustment_func = ADJUSTMENT_FUNCTIONS.get(
            variant_type, adjust_default_position
        )
        adjusted_start, adjusted_end = adjustment_func(start_position, end_position)

        if not is_supported_variant_type(variant_type, logger):
            logger.error(f"Unsupported variant type: {variant_type}")
            return None

        return (adjusted_start, adjusted_end)
    except Exception as e:
        logger.warning(f"Error adjusting positions: {e}")
        return None


# --- Counting Alternate Alleles ---


def count_alternate_alleles(variant_type, reference_allele, alternate_allele, sequence):
    """
    Counts occurrences of the alternate allele in the given sequence string from mpileup.
    """
    allele_count = 0
    allele_length = (
        len(reference_allele) if variant_type == "DEL" else len(alternate_allele)
    )

    if "NP" in variant_type:
        # MNP or DNP/TNP: Count occurrences of alternate allele as a substring
        allele_count = sequence.count(alternate_allele)
    elif variant_type == "DEL":
        # Count occurrences of deletions: represented as "-<size>"
        pattern = f"-{allele_length}"
        allele_count = len(re.findall(pattern, sequence))
    elif variant_type == "INS":
        # Count occurrences of insertions: represented as "+<size>"
        pattern = f"\\+{allele_length}"
        allele_count = len(re.findall(pattern, sequence))

    return allele_count


# --- Header Processing ---


def get_headers(maf_file, logger):
    """
    Reads the header line from a MAF file and returns a dictionary mapping
    column names to their respective indices.
    Ensures that all mandatory columns are present.
    """
    logger.debug(f"Reading headers from MAF file: {maf_file}")
    try:
        with open(maf_file, "r") as file:
            header_line = file.readline().strip()
    except Exception as e:
        logger.error(f"Failed to read MAF file {maf_file}: {e}")
        sys.exit(f"Error: Failed to read MAF file {maf_file}: {e}\n")

    header_columns = header_line.split("\t")
    header_index = {column: index for index, column in enumerate(header_columns)}
    logger.debug(f"Header columns: {header_index}")

    # Define mandatory columns
    mandatory_columns = [
        "Hugo_Symbol",
        "Gene",
        "HGVSp_Short",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Sample_Barcode",
        "Tumor_Seq_Allele2",
    ]

    missing_columns = [col for col in mandatory_columns if col not in header_index]
    if missing_columns:
        logger.error(
            f"Missing mandatory columns in MAF file {maf_file}: {', '.join(missing_columns)}"
        )
        sys.exit(
            f"Error: Missing mandatory columns in MAF file {maf_file}: {', '.join(missing_columns)}\n"
        )

    return header_index


# --- Field Extraction ---


def extract_fields(columns, header_index, logger):
    """
    Extracts necessary fields from columns based on header_index.
    Returns a dictionary of fields or None if extraction fails.
    """
    try:
        fields = {
            "chromosome": columns[header_index["Chromosome"]],
            "start_position": int(columns[header_index["Start_Position"]]),
            "reference_allele": columns[header_index["Reference_Allele"]],
            "alternate_allele": columns[header_index["Tumor_Seq_Allele2"]],
            "sample_barcode": columns[header_index["Tumor_Sample_Barcode"]],
            "variant_type": columns[header_index["Variant_Type"]],
            "end_position": int(columns[header_index["End_Position"]]),
            "Hugo_Symbol": columns[header_index["Hugo_Symbol"]],
            "Gene": columns[header_index["Gene"]],
            "HGVSp_Short": columns[header_index["HGVSp_Short"]],
        }
        return fields
    except (KeyError, IndexError, ValueError) as e:
        logger.warning(f"Error extracting fields: {e}")
        return None


# --- Variant Information Extraction ---


def construct_variant_key(fields, adjusted_positions, logger):
    """
    Constructs the variant key using chromosome, adjusted start, adjusted end, reference allele, and alternate allele.
    """
    chromosome = fields["chromosome"]
    start_position, end_position = adjusted_positions
    reference_allele = fields["reference_allele"]
    alternate_allele = fields["alternate_allele"]
    variant_key = f"{chromosome}:{start_position}:{end_position}:{reference_allele}:{alternate_allele}"
    return variant_key


def construct_variant_info(fields, logger):
    """
    Constructs the variant_info dictionary with necessary details.
    """
    variant_info = {
        "sample": fields["sample_barcode"],
        "variant_type": fields["variant_type"],
        "maf_start": fields["start_position"],
        "maf_end": fields["end_position"],
        "Hugo_Symbol": fields["Hugo_Symbol"],
        "Gene": fields["Gene"],
        "HGVSp_Short": fields["HGVSp_Short"],
    }
    return variant_info


def extract_variant_info_from_columns(columns, header_index, logger):
    """
    Extracts variant information from MAF file columns.
    Returns variant_key and variant_info if valid, else None.
    """
    fields = extract_fields(columns, header_index, logger)
    if not fields:
        return None
    adjusted_positions = adjust_positions(
        fields["variant_type"], fields["start_position"], fields["end_position"], logger
    )
    if not adjusted_positions:
        return None
    variant_key = construct_variant_key(fields, adjusted_positions, logger)
    variant_info = construct_variant_info(fields, logger)
    logger.debug(f"Extracted variant info: {variant_info}")
    return variant_key, variant_info


# --- Aggregating Variants from MAF Files ---


def add_variant_to_aggregated(variant_key, info, aggregated_variants, logger):
    """
    Adds a variant to the aggregated_variants dictionary.
    """
    if variant_key not in aggregated_variants:
        aggregated_variants[variant_key] = {
            "samples": [],
            "variant_type": info["variant_type"],
            "maf_start": info["maf_start"],
            "maf_end": info["maf_end"],
            "Hugo_Symbol": info["Hugo_Symbol"],
            "Gene": info["Gene"],
            "HGVSp_Short": info["HGVSp_Short"],
        }
    aggregated_variants[variant_key]["samples"].append(info["sample"])
    logger.debug(f"Aggregated variant {variant_key} with sample {info['sample']}.")


# --- MAF File Processing ---


def should_skip_line(line, maf_file, line_number, logger):
    """
    Determines if the current line should be skipped.
    Skips header lines and empty lines.
    """
    if is_header_line(line):
        logger.debug(f"Skipping header line in {maf_file} at line {line_number}.")
        return True
    if not is_valid_line(line, line_number, maf_file, logger):
        logger.debug(f"Skipping invalid or empty line {line_number} in {maf_file}.")
        return True
    return False


def process_line(
    line, header_index, maf_file, line_number, aggregated_variants, logger
):
    """
    Processes a single line from the MAF file.
    Extracts variant information and adds it to the aggregated_variants.
    """
    columns = split_line(line, logger, maf_file, line_number)
    if columns is None:
        return

    variant_info = extract_variant_info_from_columns(columns, header_index, logger)
    if not variant_info:
        return

    variant_key, info = variant_info
    add_variant_to_aggregated(variant_key, info, aggregated_variants, logger)


def process_maf_file(maf_file, header_index, aggregated_variants, logger):
    """
    Processes a single MAF file and updates the aggregated_variants dictionary.
    """
    logger.info(f"Processing MAF file: {maf_file}")
    try:
        with open(maf_file, "r") as file:
            for line_number, line in enumerate(file, 1):
                if should_skip_line(line, maf_file, line_number, logger):
                    continue
                process_line(
                    line,
                    header_index,
                    maf_file,
                    line_number,
                    aggregated_variants,
                    logger,
                )
    except Exception as e:
        logger.error(f"Failed to process MAF file {maf_file}: {e}")


def is_header_line(line):
    """
    Determines if a line is a header line based on the presence of 'Hugo_Symbol'.
    """
    return "Hugo_Symbol" in line


def is_valid_line(line, line_number, maf_file, logger):
    """
    Validates if a line is non-empty and not a header.
    """
    if not line.strip():
        logger.debug(f"Skipping empty line {line_number} in {maf_file}.")
        return False
    return True


def split_line(line, logger, maf_file, line_number):
    """
    Splits a line into columns. Returns None if the line is malformed.
    """
    try:
        return line.strip().split("\t")
    except Exception as e:
        logger.warning(f"Failed to split line {line_number} in {maf_file}: {e}")
        return None


def aggregate_variants_from_maf_files(maf_files, logger):
    """
    Aggregates variant information from all MAF files.
    Returns a dictionary with variant keys mapping to their aggregated information.
    """
    logger.info("Aggregating variants from MAF files.")
    aggregated_variants = {}
    for maf_file in maf_files:
        header_index = get_headers(maf_file, logger)
        process_maf_file(maf_file, header_index, aggregated_variants, logger)
    logger.info(f"Total unique variants aggregated: {len(aggregated_variants)}")
    return aggregated_variants


# --- BAM File Processing ---


def get_sample_name_from_bam(bam_file, logger):
    """
    Retrieves the sample name from a single BAM file using samtools.
    Returns the sample name or an empty string if not found.
    """
    cmd = ["samtools", "samples", bam_file]
    logger.debug(f"Running command to get sample name: {' '.join(cmd)}")
    try:
        output = (
            subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
            .decode()
            .strip()
            .split("\n")
        )
        if output:
            sample_name = output[0].split("\t")[0]
            logger.debug(
                f"Extracted sample name '{sample_name}' from BAM file '{bam_file}'."
            )
            return sample_name
        else:
            logger.warning(f"No sample name found in BAM file: {bam_file}")
            return ""
    except subprocess.CalledProcessError:
        logger.error(f"samtools command failed for BAM file: {bam_file}")
        return ""


def extract_sample_names(bam_files, logger):
    """
    Extracts sample names from BAM files using samtools.
    Returns a list of sample names.
    """
    logger.info("Extracting sample names from BAM files.")
    sample_names = []
    for bam_file in bam_files:
        sample_name = get_sample_name_from_bam(bam_file, logger)
        if not sample_name:
            logger.warning(f"Can't find sample name from BAM file {bam_file}.")
        else:
            logger.info(f"Sample '{sample_name}' inferred from BAM file '{bam_file}'.")
        sample_names.append(sample_name)
    logger.debug(f"Sample names extracted: {sample_names}")
    return sample_names


# --- Mpileup Command Construction and Execution ---


def build_mpileup_command(
    chromosome, start, end, map_quality, base_quality, bam_list_file, reference_fasta
):
    """
    Constructs the samtools mpileup command for a given variant.
    """
    return [
        "samtools",
        "mpileup",
        "-q",
        str(map_quality),
        "-Q",
        str(base_quality),
        "-d",
        "1000",
        "-b",
        bam_list_file,
        "-A",
        "-r",
        f"{chromosome}:{start}-{end}",
        "--ignore-overlaps",
        "--no-BAQ",
        "--no-output-ins",
        "--no-output-del",
        "--no-output-ends",
        "--reference",
        reference_fasta,
    ]


def run_mpileup_command(command, logger):
    """
    Executes the samtools mpileup command and captures the output.
    Returns the output as a list of lines.
    """
    logger.debug(f"Executing mpileup command: {' '.join(command)}")
    try:
        output = (
            subprocess.check_output(command, stderr=subprocess.DEVNULL)
            .decode()
            .strip()
            .split("\n")
        )
        logger.debug(f"Mpileup output received with {len(output)} lines.")
        return output
    except subprocess.CalledProcessError:
        logger.error(f"samtools mpileup command failed: {' '.join(command)}")
        return []


# --- Sample Variant Processing ---


def extract_alleles_from_key(variant_key, logger):
    """
    Extracts chromosome, start, end, reference allele, and alternate allele from the variant key.
    Returns a tuple of (chromosome, start, end, reference_allele, alternate_allele).
    """
    try:
        chromosome, start, end, reference_allele, alternate_allele = variant_key.split(
            ":"
        )
        return chromosome, int(start), int(end), reference_allele, alternate_allele
    except ValueError:
        logger.error(f"Invalid variant key format: {variant_key}")
        return None, None, None, None, None


def find_sample_index(sample_name, sample_names, logger):
    """
    Finds the index of the sample in the sample_names list.
    Returns the index or None if not found.
    """
    try:
        sample_index = sample_names.index(sample_name)
        return sample_index
    except ValueError:
        logger.error(f"Sample '{sample_name}' not found in sample list.")
        return None


def calculate_start_column(sample_index):
    """
    Calculates the starting column index for the sample in mpileup results.
    """
    return 3 + 3 * sample_index


def extract_coverage_sequence(
    mpileup_line, start_column, line_number, variant_key, sample_name, logger
):
    """
    Extracts coverage and sequence from a mpileup line.
    Returns coverage and sequence.
    """
    mpileup_columns = mpileup_line.split("\t")
    try:
        coverage = int(mpileup_columns[start_column])
        sequence = mpileup_columns[start_column + 1].upper()
        return coverage, sequence
    except (IndexError, ValueError):
        logger.warning(
            f"Invalid mpileup format for variant {variant_key} in sample {sample_name} at line {line_number + 1}."
        )
        return 0, ""


def determine_current_alt_allele(mnp_alleles, line_number, alternate_allele):
    """
    Determines the current alternate allele based on MNP alleles.
    """
    if mnp_alleles and line_number < len(mnp_alleles):
        return mnp_alleles[line_number]
    else:
        return alternate_allele


def calculate_alt_percentage(alt_count, coverage):
    """
    Calculates the alternate allele percentage.
    """
    return 100.0 * alt_count / coverage if coverage > 0 else 0


def compute_mean(values):
    """
    Computes the mean of a list of values.
    Returns the mean or 0 if the list is empty.
    """
    return statistics.mean(values) if values else 0


def check_presence_in_maf(variant_key, sample_name, aggregated_variants, logger):
    """
    Checks if the sample is present in the MAF for the given variant.
    Returns "PRESENT_IN_MAF" or "ABSENT_FROM_MAF".
    """
    samples_in_maf = aggregated_variants[variant_key]["samples"]
    presence = "PRESENT_IN_MAF" if sample_name in samples_in_maf else "ABSENT_FROM_MAF"
    return presence


def count_found_alt(alt_count_list):
    """
    Counts the number of alternate counts greater than zero.
    """
    return sum(1 for count in alt_count_list if count > 0)


def calculate_status(
    presence_in_maf, mean_alt_percentage, found_alt, mnp_alleles, logger
):
    """
    Determines the status of the variant based on presence in MAF and alternative allele percentage.
    Returns the status as a string.
    """
    status = ""
    if mnp_alleles and 0 < found_alt < len(mnp_alleles):
        status = "PART_MATCH"

    if not status:
        if presence_in_maf == "PRESENT_IN_MAF":
            status = "TRUE_POSITIVE" if mean_alt_percentage > 0 else "FALSE_POSITIVE"
        else:
            status = "TRUE_NEGATIVE" if mean_alt_percentage == 0 else "FALSE_NEGATIVE"

    logger.debug(f"Determined status: {status}")
    return status


def format_alt_percentage(mean_alt_percentage, logger):
    """
    Formats the Alt_Perc to two decimal places.
    Returns the formatted string.
    """
    if mean_alt_percentage == 0:
        formatted_alt_percentage = "0"
    else:
        formatted_alt_percentage = f"{mean_alt_percentage:.2f}"
    logger.debug(f"Formatted Alt_Perc: {formatted_alt_percentage}")
    return formatted_alt_percentage


def initialize_lists():
    """
    Initializes lists for coverage, alternate counts, and alternate percentages.
    Returns three empty lists.
    """
    return [], [], []


def process_sample_variant(
    variant_key,
    aggregated_variants,
    mpileup_results,
    sample_name,
    sample_names,
    mnp_alleles,
    variant_type,
    logger,
):
    """
    Processes the mpileup results for a single sample and determines the variant status.
    Returns a tuple of (mean_alt_count, mean_coverage, mean_alt_percentage, status).
    """
    coverage_list, alt_count_list, alt_percentage_list = initialize_lists()

    (
        chromosome,
        start,
        end,
        reference_allele,
        alternate_allele,
    ) = extract_alleles_from_key(variant_key, logger)
    if chromosome is None:
        return 0, 0, 0, "ERROR"

    sample_index = find_sample_index(sample_name, sample_names, logger)
    if sample_index is None:
        return 0, 0, 0, "ERROR"

    start_column = calculate_start_column(sample_index)

    for line_number, mpileup_line in enumerate(mpileup_results):
        coverage, sequence = extract_coverage_sequence(
            mpileup_line, start_column, line_number, variant_key, sample_name, logger
        )
        coverage_list.append(coverage)

        current_alt_allele = determine_current_alt_allele(
            mnp_alleles, line_number, alternate_allele
        )
        alt_count = count_alternate_alleles(
            variant_type, reference_allele, current_alt_allele, sequence
        )
        alt_count_list.append(alt_count)

        alt_percentage = calculate_alt_percentage(alt_count, coverage)
        alt_percentage_list.append(alt_percentage)

        logger.debug(
            f"Variant {variant_key}, Sample {sample_name}, Line {line_number + 1}: Coverage={coverage}, Alt_Count={alt_count}, Alt_Perc={alt_percentage:.2f}%"
        )

    mean_alt_count = compute_mean(alt_count_list)
    mean_coverage = compute_mean(coverage_list)
    mean_alt_percentage = compute_mean(alt_percentage_list)

    presence_in_maf = check_presence_in_maf(
        variant_key, sample_name, aggregated_variants, logger
    )
    found_alt = count_found_alt(alt_count_list)
    status = calculate_status(
        presence_in_maf, mean_alt_percentage, found_alt, mnp_alleles, logger
    )

    logger.debug(
        f"Variant {variant_key}, Sample {sample_name}: Mean_Alt_Count={mean_alt_count}, Mean_Coverage={mean_coverage}, Mean_Alt_Perc={mean_alt_percentage:.2f}%, Status={status}"
    )

    return mean_alt_count, mean_coverage, mean_alt_percentage, status


# --- Variant Processing and Output ---


def process_variants(
    aggregated_variants,
    sample_names,
    bam_list_file,
    reference_fasta,
    map_quality,
    base_quality,
    logger,
):
    """
    Iterates through each variant, runs samtools mpileup, and processes the results.
    Outputs the final results with appropriate headers.
    """
    # Define output column headers
    output_headers = [
        "Hugo_Symbol",
        "Gene",
        "HGVSp_Short",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Variant_Type",
        "Reference_Allele",
        "Tumour_Seq_Allele2",
        "Tumor_Sample_Barcode",
        "Alt_Count",
        "Tot_Count",
        "Alt_Perc",
        "Var_in_MAF",
        "Status",
    ]
    print("\t".join(output_headers))
    logger.info("Starting variant processing and output generation.")

    for variant_key in sorted(aggregated_variants.keys()):
        (
            chromosome,
            start,
            end,
            reference_allele,
            alternate_allele,
        ) = extract_alleles_from_key(variant_key, logger)
        if chromosome is None:
            continue
        start = int(start)
        end = int(end)

        logger.info(
            f"Processing variant: {chromosome}:{start}-{end} {reference_allele}->{alternate_allele}"
        )

        # Build and run the samtools mpileup command
        mpileup_command = build_mpileup_command(
            chromosome,
            start,
            end,
            map_quality,
            base_quality,
            bam_list_file,
            reference_fasta,
        )
        mpileup_output = run_mpileup_command(mpileup_command, logger)

        if not mpileup_output or not any(mpileup_output):
            logger.warning(
                f"No mpileup data for variant {chromosome}:{start}-{end} {reference_allele}->{alternate_allele}."
            )
            print(
                f"{chromosome}\t{start}\t{end}\t{reference_allele}\t{alternate_allele}\tNO_PILEUP",
                file=sys.stderr,
            )
            continue

        variant_type = aggregated_variants[variant_key]["variant_type"]
        mnp_alleles = (
            list(alternate_allele)
            if "NP" in variant_type and variant_type != "SNP"
            else []
        )

        for sample_name in sample_names:
            (
                mean_alt_count,
                mean_coverage,
                mean_alt_percentage,
                status,
            ) = process_sample_variant(
                variant_key,
                aggregated_variants,
                mpileup_output,
                sample_name,
                sample_names,
                mnp_alleles,
                variant_type,
                logger,
            )
            # Retrieve variant information
            variant_info = aggregated_variants[variant_key]
            presence_in_maf = (
                "PRESENT_IN_MAF"
                if sample_name in variant_info["samples"]
                else "ABSENT_FROM_MAF"
            )

            # Conditional formatting for Alt_Perc with fixed two decimal places
            formatted_alt_percentage = format_alt_percentage(
                mean_alt_percentage, logger
            )

            # Prepare the output line
            output_line = [
                variant_info["Hugo_Symbol"],
                variant_info["Gene"],
                variant_info["HGVSp_Short"],
                chromosome,
                variant_info["maf_start"],
                variant_info["maf_end"],
                variant_type,
                reference_allele,
                alternate_allele,
                sample_name,
                mean_alt_count,
                mean_coverage,
                formatted_alt_percentage,
                presence_in_maf,
                status,
            ]
            print("\t".join(map(str, output_line)))
            logger.debug(f"Output for sample '{sample_name}': {output_line}")


# --- Main Function ---


def main():
    """
    Main function to parse command-line arguments and orchestrate variant processing.
    """
    # Setup logging
    logger = setup_logging()
    logger.info("Variant processing script started.")

    # Parse command-line arguments
    args = parse_arguments()

    bam_list_file = args.bam_list_file
    maf_list_file = args.maf_list_file
    reference_fasta = args.reference_fasta
    map_quality = args.map_quality
    base_quality = args.base_quality

    logger.debug(
        f"Command-line arguments received: bam_list_file={bam_list_file}, maf_list_file={maf_list_file}, reference_fasta={reference_fasta}, map_quality={map_quality}, base_quality={base_quality}"
    )

    # Check required files
    check_required_files([bam_list_file, maf_list_file, reference_fasta], logger)

    # Read MAF and BAM files lists
    maf_files = read_file_list(maf_list_file, logger, file_type="MAF")
    bam_files = read_file_list(bam_list_file, logger, file_type="BAM")

    # Aggregate variants from MAF files
    aggregated_variants = aggregate_variants_from_maf_files(maf_files, logger)

    if not aggregated_variants:
        logger.error("No variants found in MAF files.")
        sys.exit("Error: No variants found in MAF files.\n")

    # Extract sample names from BAM files
    sample_names = extract_sample_names(bam_files, logger)

    if not sample_names:
        logger.error("No valid samples extracted from BAM files.")
        sys.exit("Error: No valid samples extracted from BAM files.\n")

    # Process each variant and output results
    process_variants(
        aggregated_variants,
        sample_names,
        bam_list_file,
        reference_fasta,
        map_quality,
        base_quality,
        logger,
    )

    logger.info("Variant processing script completed successfully.")


# --- Entry Point ---

if __name__ == "__main__":
    main()
