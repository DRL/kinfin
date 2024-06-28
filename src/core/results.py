import time

from core.datastore import DataFactory
from core.input import InputData
from core.utils import logger


def analyse(input_data: InputData) -> None:
    """
    Performs KinFin analysis based on the provided input data using DataFactory.

    Args:
        input_data (InputData): An instance of InputData containing input parameters and data.

    Returns:
        None

    Raises:
        Any exceptions raised by DataFactory methods.
    """
    overall_start = time.time()
    dataFactory = DataFactory(input_data)
    dataFactory.setup_dirs()
    dataFactory.analyse_clusters()
    dataFactory.aloCollection.write_tree(
        dataFactory.dirs,
        dataFactory.inputData.plot_tree,
        dataFactory.inputData.plot_format,
        dataFactory.inputData.fontsize,
    )
    dataFactory.aloCollection.compute_rarefaction_data(
        repetitions=dataFactory.inputData.repetitions,
        dirs=dataFactory.dirs,
        plotsize=dataFactory.inputData.plotsize,
        plot_format=dataFactory.inputData.plot_format,
        fontsize=dataFactory.inputData.fontsize,
    )
    dataFactory.write_output()
    overall_end = time.time()
    overall_elapsed = overall_end - overall_start
    logger.info("[STATUS] - Took %ss to run kinfin." % (overall_elapsed))
