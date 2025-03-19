import argparse
import sys

import fur_hotspot_mutations
from fur_hotspot_mutations import (
    extract_hotspot_mutations,
    mpileup_variant_filter,
    update_maf_variants,
)
from utils import constants


def main():
    # Top-level parser
    parser = argparse.ArgumentParser(
        prog=constants.PROGRAM_NAME,
        description=constants.DESCRIPTION__PROGRAM,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {fur_hotspot_mutations.__version__}",
    )

    # Command initialization
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)

    # Parsers for each command are set up by their respective modules
    _ = extract_hotspot_mutations.get_argparser(subparser=subparsers)
    _ = mpileup_variant_filter.get_argparser(subparser=subparsers)
    _ = update_maf_variants.get_argparser(subparser=subparsers)

    # Parse the arguments
    args = parser.parse_args()

    # Call the appropriate function based on the command from their respective
    # modules
    match args.command:
        case extract_hotspot_mutations.COMMAND_NAME:
            exit_code = extract_hotspot_mutations.main(args)
        case mpileup_variant_filter.COMMAND_NAME:
            exit_code = mpileup_variant_filter.main(args)
        case update_maf_variants.COMMAND_NAME:
            exit_code = update_maf_variants.main(args)
        case _:
            # Unlikley to see this error as the subparsers as
            # parser.parse_args() will catch unknown commands
            raise argparse.ArgumentError(
                argument=None,
                message=f"Invalid command - {args.command}",
            )

    # Exit with the exit code from the called function
    if exit_code is not None:
        sys.exit(exit_code)

    return None
