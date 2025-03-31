import typing as t
import logging
import sys

from utils import constants


def setup_logging(level: str = "INFO", name: t.Optional[str] = None) -> logging.Logger:
    """
    Configures the logging settings.
    Logs are output to stderr with a specific format and levels.
    """
    if name is None:
        name = constants.LOGGER_NAME

    logging_level = _get_logging_level_enum(level)

    logger = logging.getLogger(name)
    logger.setLevel(logging_level)

    # Create handlers
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)  # Capture all levels in the handler

    # Create formatter and add it to the handler
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)

    # Add handler to the logger
    logger.addHandler(handler)

    return logger


def update_logger_level(logger: logging.Logger, level: str) -> None:
    """
    Update the logger level.
    """
    logging_level = _get_logging_level_enum(level)

    logger.setLevel(logging_level)
    return None


def _get_logging_level_enum(level: str) -> int:
    """
    Get the logging level enum from the logging module.
    """
    logging_level = getattr(logging, level.strip().upper())
    if not isinstance(logging_level, int):
        raise ValueError(f"Invalid log level: {level}")
    return logging_level
