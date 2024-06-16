import gzip
import os
import sys
from typing import Any, Generator


def progress(iteration: int, steps: int | float, max_value: int) -> None:
    if int(iteration) == int(max_value):
        sys.stdout.write("\r")
        print("[PROGRESS] \t- %d%%" % (100))
    elif int(iteration) % int(steps + 1) == 0:
        sys.stdout.write("\r")
        print(
            "[PROGRESS] \t- %d%%" % (float(int(iteration) / int(max_value)) * 100),
            end=" ",
        )
        sys.stdout.flush()
    else:
        pass


def check_file(filepath: str | None, install_kinfin: bool = False) -> None:
    """
    Check if a file exists.

    Args:
        filepath (str): Path to the file to be checked.

    Raises:
        FileNotFoundError: If the file does not exist.
    """

    if filepath is not None and not os.path.isfile(filepath):
        error_msg = f"[ERROR] - file {filepath} not found."
        if install_kinfin:
            error_msg += " Please run the install script to download kinfin."
        raise FileNotFoundError(error_msg)


def yield_file_lines(filepath: str) -> Generator[str, Any, None]:
    """
    Args:
        filepath (str): Path to the file.

    Yields:
        str: Each line from the file.
    """
    check_file(filepath)
    if filepath.endswith(".gz"):
        with gzip.open(filepath, "rb") as fh:
            for line in fh:
                line = line.decode("utf-8")
                if line.startswith("nodesDB.txt"):
                    line = "#%s" % line.split("#")[1]
                yield line.rstrip("\n")
    else:
        with open(filepath) as fh:
            for line in fh:
                yield line.rstrip("\n")
