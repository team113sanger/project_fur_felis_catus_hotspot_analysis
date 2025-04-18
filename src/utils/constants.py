"""
This module contains string constants used in the CLI, centralizing
them in one place for easy maintenance.
"""
# Logger name
LOGGER_NAME = "FurHotspotMutations"

# Program name
#
# The name of the program that appears in the unified CLI entrypoint help message.
PROGRAM_NAME: str = "fur_hotspot_mutations"


# Command names
#
# Command names appear in the help message and are used to determine which
# subcommand to run.
COMMAND_NAME__EXTRACT_HOTSPOT_MUTATIONS: str = "extract_hotspot_mutations"
COMMAND_NAME__MPILEUP_VARIANT_FILTER: str = "mpileup_variant_filter"
COMMAND_NAME__UPDATE_MAF_VARIANTS: str = "update_maf_variants"
COMMAND_NAME__CHECK_MPILEUP_SAMPLES_FROM_MAF: str = "check_mpileup_samples_from_maf"

# Program & Command descriptions
#
# Text that appears in the help message that describes a command or the program
# itself.
DESCRIPTION__PROGRAM: str = (
    "A toolkit for identifying hotspot mutations in tumor cohorts "
    "and rescuing false-negative variant calls from Caveman and Pindel. "
    "The program processes MAF files to evaluate mutations against specific criteria, "
    "rescue false-negatives, remove false-positives, and reclassify suspected germline mutations."
)
DESCRIPTION__EXTRACT_HOTSPOT_MUTATIONS: str = (
    "Identify and extract hotspot mutations from a Mutation Annotation Format (MAF) file. "
    "Hotspot mutations are defined as identical mutations occurring at the same site in more than N samples. "
    "The extracted hotspot mutations are saved to a new MAF file."
)
DESCRIPTION__MPILEUP_VARIANT_FILTER: str = (
    "Uses mpileup data to identify questionable variant calls. "
    "Uses thresholds determined by the user to decide whether a given mutation "
    "should be included/excluded in the final callset."
)
DESCRIPTION__UPDATE_MAF_VARIANTS: str = (
    "This command takes a MAF file as input, adds/removes variants "
    "according to a variant file, and then writes the updated variants "
    "to a new MAF file."
)
DESCRIPTION__CHECK_MPILEUP_SAMPLES_FROM_MAF: str = (
    "Process MAF and BAM files to analyze genetic variants using samtools mpileup."
)


# Command short-help
#
# Text that appears in the program help that concisely describes the command.
SHORT_HELP__EXTRACT_HOTSPOT_MUTATIONS: str = (
    "Identify and extract hotspot mutations from a MAF file."
)
SHORT_HELP__MPILEUP_VARIANT_FILTER: str = "Filter variant calls using mpileup data."
SHORT_HELP__UPDATE_MAF_VARIANTS: str = "Add or remove variants in a MAF file."
SHORT_HELP__CHECK_MPILEUP_SAMPLES_FROM_MAF: str = (
    "Process MAF and BAM files using samtools mpileup."
)
