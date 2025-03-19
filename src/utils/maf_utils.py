import csv
from pathlib import Path
from typing import Generator, List, Optional

from utils.logging_utils import setup_logging

LOGGER = setup_logging()

# CONSTANTS
EXPECTED_COLUMNS = {
    "hugo_symbol",
    "center",
    "ncbi_build",
    "chromosome",
    "start_position",
    "end_position",
    "variant_classification",
    "variant_type",
    "reference_allele",
    "tumor_seq_allele2",
}


# HELPERS


# File Reading Helpers
def _read_file(file_path: Path) -> Generator[List[str], None, None]:
    """
    Reads a tab-delimited file line by line and yields rows as lists of strings.

    Ignores comment lines (starting with '#') and empty lines.

    Args:
        file_path (Path): The path to the file to read.

    Yields:
        List[str]: A list of strings representing the columns in a non-comment row.

    Raises:
        Exception: If the file cannot be read or accessed.
    """
    try:
        with file_path.open("r", newline="") as file:
            reader = csv.reader(file, delimiter="\t")
            for row in reader:
                if row and not row[0].startswith("#"):  # Skip empty and comment lines
                    yield row
    except Exception as e:
        LOGGER.exception(f"Error reading {file_path}: {e}")
        raise


# Validation Helpers
def _validate_header(header: List[str], file_path: Path) -> Optional[List[str]]:
    """
    Validates the header row of a file to ensure it contains all required columns.

    Normalizes column names by stripping whitespace and converting to lowercase.

    Args:
        header (List[str]): A list of column names from the file header.
        file_path (Path): The path to the file being validated.

    Returns:
        Optional[List[str]]: The normalized header if valid; None otherwise.

    Logs:
        A warning if required columns are missing.
    """
    normalized_header = [col.strip().lower() for col in header]
    missing_columns = EXPECTED_COLUMNS - set(normalized_header)
    if missing_columns:
        LOGGER.warning(
            f"{file_path} is missing the following expected MAF columns: {', '.join(missing_columns)}"
        )
        return None
    return normalized_header


def _validate_positions(
    row: List[str], line_number: int, file_path: Path, header_index: dict
) -> bool:
    """
    Validates that the 'start_position' and 'end_position' columns contain valid integers
    and that 'start_position' is less than or equal to 'end_position'.

    Args:
        row (List[str]): A list of column values for the current row.
        line_number (int): The line number in the file (1-based).
        file_path (Path): The path to the file being validated.
        header_index (dict): A mapping of column names to their indices.

    Returns:
        bool: True if the positions are valid; False otherwise.

    Logs:
        A warning if positions are invalid or non-integer.
    """
    try:
        start_pos = int(row[header_index["start_position"]])
        end_pos = int(row[header_index["end_position"]])
        if start_pos > end_pos:
            LOGGER.warning(
                f"Invalid positions on line {line_number} in {file_path}: "
                f"Start_Position ({start_pos}) > End_Position ({end_pos})."
            )
            return False
    except ValueError:
        LOGGER.warning(
            f"Non-integer positions found on line {line_number} in {file_path}: "
            f"Start_Position or End_Position is not an integer."
        )
        return False
    return True


def _validate_data_lines(
    lines: Generator[List[str], None, None], header: List[str], file_path: Path
) -> bool:
    """
    Validates all data rows in the file.

    Ensures the start and end positions are valid for each row.

    Args:
        lines (Generator[List[str], None, None]): A generator yielding rows from the file.
        header (List[str]): The normalized header of the file.
        file_path (Path): The path to the file being validated.

    Returns:
        bool: True if all rows are valid; False otherwise.
    """
    LOGGER.info("Validating MAF data lines ...")
    header_index = {col: idx for idx, col in enumerate(header)}
    for line_number, row in enumerate(lines, start=2):  # Start after header
        if not _validate_positions(row, line_number, file_path, header_index):
            return False
    LOGGER.info("Successfully validated MAF data lines.")
    return True


# Duplicate Line Checking Helpers
def _check_for_duplicates(file_path: Path) -> bool:
    """
    Checks the file for duplicate rows based on their content.

    Args:
        file_path (Path): The path to the file to check.

    Returns:
        bool: True if no duplicates are found; False otherwise.

    Logs:
        A warning if duplicates are detected or if the file has no valid lines.
    """
    LOGGER.info(f"Checking {file_path} for duplicate lines ...")
    result = _detect_duplicates(file_path)
    if result is None:
        LOGGER.warning(
            f"The file {file_path} contains no valid lines (empty or only comments)."
        )
        return False
    if result:
        LOGGER.warning(f"The file {file_path} contains duplicate lines.")
        return False
    LOGGER.info(f"No duplicates found in {file_path}.")
    return True


def _detect_duplicates(file_path: Path) -> Optional[bool]:
    """
    Detects duplicate rows in the file.

    Args:
        file_path (Path): The path to the file to check.

    Returns:
        Optional[bool]: True if duplicates are found, False if no duplicates,
                        or None if the file has no valid lines.
    """
    seen = {}
    lines_checked = 0

    for line_number, row in enumerate(_read_file(file_path), start=1):
        lines_checked += 1
        if _is_duplicate(row, seen, line_number, file_path):
            return True
    return _handle_no_lines_checked(lines_checked, file_path)


def _is_duplicate(
    row: List[str], seen: dict, line_number: int, file_path: Path
) -> bool:
    """
    Checks if a specific row is a duplicate of a previously seen row.

    Args:
        row (List[str]): A list of column values for the current row.
        seen (dict): A dictionary mapping row content to its first occurrence line number.
        line_number (int): The current line number in the file (1-based).
        file_path (Path): The path to the file being checked.

    Returns:
        bool: True if the row is a duplicate; False otherwise.

    Logs:
        A warning if a duplicate is found.
    """
    line_str = "\t".join(row)
    if line_str in seen:
        LOGGER.warning(
            f"Duplicate line found in {file_path}: Line {seen[line_str]} and Line {line_number}."
        )
        return True
    seen[line_str] = line_number
    return False


def _handle_no_lines_checked(lines_checked: int, file_path: Path) -> Optional[bool]:
    """
    Handles the case where no lines were checked (empty or comment-only file).

    Args:
        lines_checked (int): The total number of lines checked.
        file_path (Path): The path to the file being checked.

    Returns:
        Optional[bool]: None if no lines were checked; False otherwise.

    Logs:
        A warning if the file contains no valid lines.
    """
    if lines_checked == 0:
        LOGGER.warning(f"The file {file_path} contains no valid lines.")
        return None
    return False


# Column Consistency Checking Helpers
def _check_column_consistency(file_path: Path) -> bool:
    """
    Checks if all rows in the file have the same number of columns.

    Args:
        file_path (Path): The path to the file to check.

    Returns:
        bool: True if all rows have consistent column counts; False otherwise.

    Logs:
        A warning if column counts are inconsistent or if the file is empty.
    """
    LOGGER.info(f"Checking {file_path} for consistent column counts ...")
    try:
        with file_path.open("r", newline="") as file:
            reader = csv.reader(file, delimiter="\t")
            rows = list(reader)
            if not rows:
                LOGGER.warning(f"The file {file_path} is empty.")
                return False
            expected_columns = len(rows[0])
            if all(len(row) == expected_columns for row in rows):
                LOGGER.info(
                    f"Successfully checked {file_path} for consistent column counts."
                )
                return True
            LOGGER.warning(f"The file {file_path} has inconsistent column counts.")
            return False
    except Exception as e:
        LOGGER.exception(f"Error checking for consistent columns in {file_path}: {e}")
        return False


# MAIN FUNCTIONS


def _check_file_structure(file_path: Path) -> bool:
    """
    Performs initial structure checks on the file:
    - Duplicate lines check.
    - Column consistency check.

    Args:
        file_path (Path): The path to the file.

    Returns:
        bool: True if the file passes the structural checks; False otherwise.
    """
    if not _check_for_duplicates(file_path):
        return False
    if not _check_column_consistency(file_path):
        return False
    return True


def _validate_header_and_data(file_path: Path) -> bool:
    """
    Validates the file's header and data rows.

    Args:
        file_path (Path): The path to the file.

    Returns:
        bool: True if header and data rows are valid; False otherwise.
    """
    lines = _read_file(file_path)
    header = next(lines, None)  # Get the first line (header)
    if not header:
        LOGGER.warning(f"No header found in {file_path}.")
        return False
    normalized_header = _validate_header(header, file_path)
    if not normalized_header:
        return False
    return _validate_data_lines(lines, normalized_header, file_path)


def is_maf_format(file_path: Path) -> bool:
    """
    Determines whether a file adheres to the Mutation Annotation Format (MAF) standard.

    This function orchestrates the following checks:
    1. File structure (duplicates, column consistency).
    2. Header and data row validation.

    Args:
        file_path (Path): The path to the file to validate.

    Returns:
        bool: True if the file is valid MAF format; False otherwise.
    """
    try:
        if not _check_file_structure(file_path):
            return False
        if not _validate_header_and_data(file_path):
            return False
        return True
    except Exception as e:
        LOGGER.exception(f"Error validating MAF format for {file_path}: {e}")
        return False
