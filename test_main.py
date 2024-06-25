import os
import pytest


def get_file_pairs():
    dir1 = "result/advanced.cli.old"
    dir2 = "result/advanced.cli.new"

    assert os.path.exists(dir1), f"Directory '{dir1}' does not exist"
    assert os.path.exists(dir2), f"Directory '{dir2}' does not exist"

    files1 = get_files(dir1)
    files2 = get_files(dir2)

    assert len(files1) == len(
        files2
    ), "Directories do not contain the same number of files"

    file_pairs = []
    for file1, file2 in zip(files1, files2):
        if file1.endswith(".txt"):
            file_pairs.append((os.path.join(dir1, file1), os.path.join(dir2, file2)))
    return file_pairs


def get_files(directory):
    """
    Recursively get all files in a directory
    """
    file_list = []
    for root, _, files in os.walk(directory):
        for file in files:
            relative_path = os.path.relpath(os.path.join(root, file), directory)
            file_list.append(relative_path)
    return file_list


def compare_files(file1, file2):
    """
    Compare files based on their types
    """
    if file1.endswith(".txt"):
        return check_is_mismatch(file1, file2)
    else:
        return False


def check_is_mismatch(file1, file2):
    """
    Compare each line of two text files
    """
    with open(file1, "r") as f1, open(file2, "r") as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

        # Remove empty lines and strip whitespace
        lines1 = [line.strip() for line in lines1 if line.strip()]
        lines2 = [line.strip() for line in lines2 if line.strip()]

        # Sort lines
        lines1.sort()
        lines2.sort()

        # Compare sorted lines
        if lines1 == lines2:
            return False
        else:
            return True


@pytest.mark.parametrize("file1, file2", get_file_pairs())
def test_compare_files(file1, file2):
    mismatch = compare_files(file1, file2)
    assert not mismatch, f"Files '{file1}' and '{file2}' have mismatches"
