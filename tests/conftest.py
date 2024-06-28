import os
from typing import List, Tuple


def pytest_addoption(parser) -> None:
    """Add argument to take path to generated and expected output directories"""
    parser.addoption(
        "--generated",
        action="store",
        help="Path to the generated output directory",
    )
    parser.addoption(
        "--expected",
        action="store",
        help="Path to the expected output directory",
    )


def pytest_generate_tests(metafunc) -> None:
    """Generates test for each file"""
    if "gen_file" in metafunc.fixturenames and "exp_file" in metafunc.fixturenames:
        file_pairs = get_file_pairs(metafunc.config)
        metafunc.parametrize("gen_file,exp_file", file_pairs)


def get_file_pairs(config) -> List[Tuple[str, str]]:
    """Get tuple of generate result file vs expected result file to compare"""
    generated = config.getoption("generated")
    expected = config.getoption("expected")

    assert os.path.exists(generated), f"Directory '{generated}' does not exist"
    assert os.path.exists(expected), f"Directory '{expected}' does not exist"

    files1: List[str] = get_files(generated)
    files2: List[str] = get_files(expected)

    assert len(files1) == len(
        files2
    ), "Directories do not contain the same number of files"

    file_pairs: List[Tuple[str, str]] = []
    for gen_file, exp_file in zip(files1, files2):
        if gen_file.endswith(".txt"):
            file_pairs.append(
                (os.path.join(generated, gen_file), os.path.join(expected, exp_file))
            )
    return file_pairs


def get_files(directory) -> List[str]:
    """
    Recursively get all files in a directory
    """
    file_list = []
    for root, _, files in os.walk(directory):
        for file in files:
            relative_path = os.path.relpath(os.path.join(root, file), directory)
            file_list.append(relative_path)
    return file_list
