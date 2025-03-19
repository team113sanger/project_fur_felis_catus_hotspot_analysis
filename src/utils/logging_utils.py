import logging


def setup_logging(level: str = "DEBUG") -> None:
    logging_level = getattr(logging, level.strip().upper())
    if not isinstance(logging_level, int):
        raise ValueError(f"Invalid log level: {level}")
    logging.basicConfig(
        level=logging_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )
    logging.info("Logging setup complete with level: %s", level)
    return None
