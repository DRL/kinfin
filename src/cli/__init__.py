from core.input import InputData
from core.results import analyse


def run_cli(args: InputData) -> None:
    """
    Run the command-line interface to perform analysis based on the provided input data.

    Args:
        args (InputData): An instance of InputData containing input parameters and data.

    Returns:
        None
    """
    analyse(args)
