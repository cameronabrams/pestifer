# Author: ChatGPT 5
import os
_MODE: bool | None = None   # None = auto (fallback to env); True/False = explicit
_EXAMPLE_ID: int | None = None  # for reporting which example is being processed

def set_generate_gold(val: bool) -> None:
    """Called by the pytest plugin to force the mode during tests."""
    global _MODE
    _MODE = bool(val)

def generate_gold() -> bool:
    """Call this from your package code to decide save vs compare."""
    if _MODE is not None:
        return _MODE
    return os.getenv("PESTIFER_PYTEST_GENERATE_GOLD") == "1"

def set_example_id(example_id: int) -> None:
    """Set the example ID being processed."""
    global _EXAMPLE_ID
    _EXAMPLE_ID = example_id

def report_example_id() -> int | None:
    """Report the example ID being processed, if any."""
    if _EXAMPLE_ID is not None:
        return _EXAMPLE_ID
    return int(os.getenv("PESTIFER_PYTEST_GENERATE_GOLD_EXAMPLE_ID", 0))