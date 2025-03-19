import argparse

import fur_hotspot_mutations
from fur_hotspot_mutations import extract_hotspot_mutations
from utils import constants, logging_utils


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
    _ = subparsers.add_parser("two", help="command two help")  # TODO: replace

    # Parse the arguments
    args = parser.parse_args()

    # Call the appropriate function based on the command from their respective
    # modules
    logging_utils.setup_logging()
    match args.command:
        case extract_hotspot_mutations.COMMAND_NAME:
            extract_hotspot_mutations.main(args)
        case "two":
            do_two(args)
        case _:
            # Unlikley to see this error as the subparsers as
            # parser.parse_args() will catch unknown commands
            raise argparse.ArgumentError(
                argument=None,
                message=f"Invalid command - {args.command}",
            )

    return None


def do_two(args: argparse.Namespace) -> None:  # TODO: remove
    print("Doing two")
    print(args)
