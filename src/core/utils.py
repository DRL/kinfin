import gzip
import logging
import os
import sys
from math import log, sqrt
from typing import Any, Generator, List, Optional, Tuple

import scipy

# setup logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(message)s", "%Y-%m-%d %H:%M:%S")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# log_file = "app.log"
# file_handler = logging.FileHandler(log_file)
# file_handler.setFormatter(formatter)
# logger.addHandler(file_handler)


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


def read_fasta_len(fasta_file: str) -> Generator[Tuple[str, int], Any, None]:
    """
    Generator function to parse a FASTA file and yield tuples of header and sequence length.

    Args:
    - fasta_file (str): Path to the FASTA file to be parsed.

    Yields:
    Tuple[str, int]: A tuple containing the header and the length of the sequence.

    Raises:
    FileNotFoundError: If the specified FASTA file does not exist.
    """
    check_file(fasta_file)
    with open(fasta_file) as fh:
        print(f"[STATUS]\t - Parsing FASTA {fasta_file}")
        header: str = ""
        seqs: List[str] = []
        for line in fh:
            if line[0] == ">":
                if header:
                    header = (
                        header.replace(":", "_")
                        .replace(",", "_")
                        .replace("(", "_")
                        .replace(")", "_")
                    )  # orthofinder replaces chars
                    yield header, len("".join(seqs))
                header, seqs = (
                    line[1:-1].split()[0],
                    [],
                )  # Header is split at first whitespace
            else:
                seqs.append(line[:-1])
        header = (
            header.replace(":", "_")
            .replace(",", "_")
            .replace("(", "_")
            .replace(")", "_")
        )  # orthofinder replaces chars
        yield header, len("".join(seqs))


def median(lst) -> float:
    list_sorted = sorted(lst)
    list_length = len(lst)
    index = (list_length - 1) // 2
    if list_length % 2:
        return list_sorted[index] / 1.0
    else:
        return (list_sorted[index] + list_sorted[index + 1]) / 2.0


def mean(lst) -> float:
    if lst:
        return float(sum(lst)) / len(lst)
    else:
        return 0.0


def sd(lst, population=True) -> float:
    n = len(lst)
    differences = [x_ - mean(lst) for x_ in lst]
    sq_differences = [d**2 for d in differences]
    ssd = sum(sq_differences)
    if population is True:
        variance = ssd / n
    else:
        variance = ssd / (n - 1)
    sd_result = sqrt(variance)
    return sd_result


def statistic(
    count_1: List[int],
    count_2: List[int],
    test: str,
    min_proteomes: int,
) -> Tuple[
    Optional[float],
    Optional[float],
    Optional[float],
    Optional[float],
]:
    pvalue: Optional[float] = None
    log2_mean: Optional[float] = None
    mean_count_1: Optional[float] = None
    mean_count_2: Optional[float] = None

    implicit_count_1: List[float] = [count for count in count_1 if count > 0]
    implicit_count_2: List[float] = [count for count in count_2 if count > 0]

    if len(implicit_count_1) < min_proteomes or len(implicit_count_2) < min_proteomes:
        return None, None, None, None

    mean_count_1 = mean(implicit_count_1)
    mean_count_2 = mean(implicit_count_2)
    log2_mean = log(mean_count_1 / mean_count_2, 2)

    if (
        len(set(implicit_count_1)) == 1
        and len(set(implicit_count_2)) == 1
        and set(implicit_count_1) == set(implicit_count_2)
    ):  # equal
        pvalue = 1.0
    elif test == "welch":
        # try:
        # Welch's t-test
        pvalue = scipy.stats.ttest_ind(
            implicit_count_1,
            implicit_count_2,
            equal_var=False,
        )[1]

        if pvalue != pvalue:  # testing for "nan"
            pvalue = 1.0
    elif test == "mannwhitneyu":
        try:
            pvalue = scipy.stats.mannwhitneyu(
                implicit_count_1,
                implicit_count_2,
                alternative="two-sided",
            )[1]
        except ValueError:  # throws ValueError when all numbers are equal
            pvalue = 1.0
    elif test == "ttest":
        # try:
        pvalue = scipy.stats.ttest_ind(implicit_count_1, implicit_count_2)[1]  # t-test
        if pvalue != pvalue:  # testing for "nan"
            pvalue = 1.0
    elif test == "ks":
        # H0 that they are drawn from the same distribution
        pvalue = scipy.stats.ks_2samp(implicit_count_1, implicit_count_2)[1]
        if pvalue != pvalue:  # testing for "nan"
            pvalue = 1.0
    elif test == "kruskal":
        # H0 is that population median is equal
        pvalue = scipy.stats.kruskal(implicit_count_1, implicit_count_2)[1]
        if pvalue != pvalue:  # testing for "nan"
            pvalue = 1.0
    else:
        pass

    return pvalue, log2_mean, mean_count_1, mean_count_2
