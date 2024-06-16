import argparse

from core.config import SUPPORTED_PLOT_FORMATS, SUPPORTED_TAXRANKS, SUPPORTED_TESTS


# TODO : --plotsize should take a tuple
# TODO : --taxranks should take multiple inputs
def parse_args() -> None:
    """Parse command-line arguments.

    Args:
        nodesdb_f (str): filepath of nodesdb_f.
        pfam_mapping_f (str): filepath of pfam_mapping_f.
        ipr_mapping_f (str): filepath of ipr_mapping_f.
        go_mapping_f (str): filepath of go_mapping_f.

    Returns:
        ServeArgs or InputData: Parsed arguments based on the command.

    Raises:
        SystemExit: If an invalid command is provided.
    """

    parser = argparse.ArgumentParser(
        description="Kinfin proteome cluster analysis tool"
    )

    subparsers = parser.add_subparsers(title="command", required=True, dest="command")
    api_parser = subparsers.add_parser("serve", help="Start the server")
    api_parser.add_argument(
        "-p",
        "--port",
        type=int,
        default=8000,
        help="Port number for the server (default: 8000)",
    )

    cli_parser = subparsers.add_parser("analyse", help="Perform analysis")

    # Required Arguments
    required_group = cli_parser.add_argument_group("Required Arguments")
    required_group.add_argument(
        "-g",
        "--cluster_file",
        help="OrthologousGroups.txt produced by OrthoFinder",
        required=True,
    )
    required_group.add_argument(
        "-c", "--config_file", help="Config file (in CSV format)", required=True
    )
    required_group.add_argument(
        "-s",
        "--sequence_ids_file",
        help="SequenceIDs.txt used in OrthoFinder",
        required=True,
    )

    # Other Files
    other_files_group = cli_parser.add_argument_group("Other Files")
    other_files_group.add_argument(
        "-p", "--species_ids_file", help="SpeciesIDs.txt used in OrthoFinder"
    )
    other_files_group.add_argument(
        "-f",
        "--functional_annotation",
        help="Mapping of ProteinIDs to GO/IPRS/SignalP/Pfam (can be generated through 'iprs_to_table.py')",
    )
    other_files_group.add_argument("-a", "--fasta_dir", help="Directory of FASTA files")
    other_files_group.add_argument(
        "-t",
        "--tree_file",
        help="Tree file in Newick format (taxon names must be the same as TAXON in config file)",
    )

    # General Options
    general_group = cli_parser.add_argument_group("General Options")
    general_group.add_argument("-o", "--output_path", help="Output prefix")
    general_group.add_argument(
        "--infer_singletons",
        help="Absence of proteins in clustering is interpreted as singleton (based on SequenceIDs.txt)",
        action="store_true",
    )
    general_group.add_argument(
        "--plot_tree",
        help="Plot PDF of annotated phylogenetic tree (requires -t, full ETE3 installation and X-server/xvfb-run)",
        action="store_true",
    )
    general_group.add_argument(
        "--min_proteomes",
        help="Required number of proteomes in a taxon-set to be used in rarefaction/representation-test computations [default: 2]",
        default=2,
        type=int,
    )
    general_group.add_argument(
        "--test",
        help="Test to be used in representation-test computations [default: mannwhitneyu]. Options: ttest, welch, mannwhitneyu, ks, kruskal",
        default="mannwhitneyu",
        choices=SUPPORTED_TESTS,
    )
    general_group.add_argument(
        "-r",
        "--taxranks",
        help="Taxonomic ranks to be inferred from TaxIDs in config file [default: phylum,order,genus]",
        # TODO : Add SUPPORTED_TAXRANKS here
        default=["phylum", "order", "genus"],
        nargs="+",
        choices=SUPPORTED_TAXRANKS,
    )
    general_group.add_argument(
        "--repetitions",
        help="Number of repetitions for rarefaction curves [default: 30]",
        default=30,
        type=int,
    )

    # Fuzzy Orthology Groups
    fuzzy_group = cli_parser.add_argument_group("Fuzzy Orthology Groups")
    fuzzy_group.add_argument(
        "-n",
        "--target_count",
        help="Target number of copies per proteome [default: 1]",
        default=1,
        type=int,
    )
    fuzzy_group.add_argument(
        "-x",
        "--target_fraction",
        help="Min proportion of proteomes at target_count [default: 0.75]",
        default=0.75,
        type=float,
    )
    fuzzy_group.add_argument(
        "--min",
        help="Min count of proteins for proteomes outside of target_fraction [default: 0]",
        default=0,
        type=int,
    )
    fuzzy_group.add_argument(
        "--max",
        help="Max count of proteins for proteomes outside of target_fraction [default: 20]",
        default=20,
        type=int,
    )

    plotting_group = cli_parser.add_argument_group("Plotting Options")
    plotting_group.add_argument(
        "--fontsize", help="Fontsize for plots [default: 18]", default=18, type=int
    )
    plotting_group.add_argument(
        "--plotsize",
        help="Size (WIDTH,HEIGHT) for plots [default: 24,12]",
        default=(24, 12),
        nargs=2,
    )
    plotting_group.add_argument(
        "--plot_format",
        help="Plot formats [default: pdf]",
        default="pdf",
        choices=SUPPORTED_PLOT_FORMATS,
    )
