def compare_files(gen_file, exp_file):
    """
    Compare files based on their types
    """
    if gen_file.endswith(".txt"):
        return check_is_mismatch(gen_file, exp_file)
    else:
        return False


def check_is_mismatch(gen_file, exp_file):
    """
    Compare each line of two text files
    """
    with open(gen_file, "r") as f1, open(exp_file, "r") as f2:
        gen_lines = f1.readlines()
        exp_lines = f2.readlines()
    # Remove empty lines and strip whitespace
    gen_lines = [line.strip() for line in gen_lines if line.strip()]
    exp_lines = [line.strip() for line in exp_lines if line.strip()]
    # Sort lines
    gen_lines.sort()
    exp_lines.sort()
    # Compare sorted lines
    return gen_lines != exp_lines


def test_compare_files(gen_file, exp_file):
    mismatch = compare_files(gen_file, exp_file)
    assert not mismatch, f"Files '{gen_file}' and '{exp_file}' have mismatches"
