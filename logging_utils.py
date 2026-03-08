import logging
import os
from typing import Optional


def get_conversation_logger(log_path: Optional[str] = None) -> logging.Logger:
    """
    Return a process-wide logger for agent conversations.

    All processes (client, orchestrator, sub-agents) should use this so that
    they append to the same log file with a consistent format.

    If log_path is None, it will use:
      - AGENT_LOG_PATH env var, or
      - 'agent_conversation.log' in the current working directory.
    """
    logger = logging.getLogger("conversation")

    # Avoid configuring the logger multiple times.
    if getattr(logger, "_a2a_configured", False):
        return logger

    if log_path is None:
        log_path = os.getenv("AGENT_LOG_PATH", "agent_conversation.log")

    handler = logging.FileHandler(log_path, encoding="utf-8")
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    )

    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    logger.propagate = False

    # Mark as configured to prevent duplicate handlers.
    setattr(logger, "_a2a_configured", True)

    return logger

