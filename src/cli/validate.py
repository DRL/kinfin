import sys

from core.utils import check_file, logger


def validate_cli_args(args) -> None:
    """Validate cli input arguments.

    This function checks if all required files exist and if the arguments meet specific conditions.

    Args:
        args (InputData): Input arguments as a named tuple.

    Raises:
        SystemExit: If there are any validation errors, exits the program with error messages.
    """

    error_msgs = []

    try:
        check_file(args.cluster_file)
    except FileNotFoundError as e:
        error_msgs.append(str(e))
    try:
        if not isinstance(args.config_file, str):
            raise ValueError("[ERROR] - Invalid config file data")

        check_file(args.config_file)
    except (FileNotFoundError, ValueError) as e:
        error_msgs.append(str(e))
    try:
        check_file(args.sequence_ids_file)
    except FileNotFoundError as e:
        error_msgs.append(str(e))
    try:
        check_file(args.species_ids_file)
    except FileNotFoundError as e:
        error_msgs.append(str(e))
    try:
        check_file(args.tree_file)
    except FileNotFoundError as e:
        error_msgs.append(str(e))
    try:
        check_file(args.functional_annotation)
    except FileNotFoundError as e:
        error_msgs.append(str(e))

    if args.fasta_dir and not args.species_ids_file:
        error_msgs.append(
            "[ERROR] : You have provided a FASTA-dir using '--fasta-dir'. Please also provide a Species-ID file using ('--species_ids_file')."
        )

    if args.target_count < 0:
        error_msgs.append(
            f"[ERROR] : --target_count {args.target_count} must be greater than 0"
        )

    if args.target_fraction < 0 or args.target_fraction > 1:
        error_msgs.append(
            f"[ERROR] : --target_fraction {args.target_fraction} is not between 0.0 and 1.0"
        )

    if args.min > args.max:
        error_msgs.append(
            f"[ERROR] : --min {args.min} is greater than --max {args.max}"
        )

    if not args.repetitions > 0:
        error_msgs.append(
            "[ERROR] : Please specify a positive integer for the number of repetitions for the rarefaction curves"
        )

    if not args.min_proteomes > 0:
        error_msgs.append(
            "[ERROR] : Please specify a positive integer for the minimum number of proteomes to consider for computations"
        )

    if error_msgs:
        logger.error("\n".join(error_msgs))
        sys.exit(1)
