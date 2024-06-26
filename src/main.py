#!/usr/bin/env python3

import os
import sys

from api import run_server
from cli import run_cli
from cli.commands import parse_args
from core.input import InputData, ServeArgs
from core.utils import check_file, logger

if __name__ == "__main__":

    # Without these files, application won't start
    base_dir = os.getcwd()
    nodesdb_f = os.path.join(base_dir, "data/nodesdb.txt")
    pfam_mapping_f = os.path.join(base_dir, "data/Pfam-A.clans.tsv.gz")
    ipr_mapping_f = os.path.join(base_dir, "data/entry.list")
    go_mapping_f = os.path.join(base_dir, "data/interpro2go")

    try:
        check_file(nodesdb_f, install_kinfin=True)
        check_file(pfam_mapping_f, install_kinfin=True)
        check_file(ipr_mapping_f, install_kinfin=True)
        check_file(go_mapping_f, install_kinfin=True)
    except FileNotFoundError as e:
        logger.error(e)
        sys.exit(1)

    args = parse_args(nodesdb_f, pfam_mapping_f, ipr_mapping_f, go_mapping_f)

    if isinstance(args, ServeArgs):
        # run the api server
        cluster_f = os.environ.get("CLUSTER_FILE_PATH")
        sequence_ids_f = os.environ.get("SEQUENCE_IDS_FILE_PATH")
        taxon_idx_mapping_file = os.environ.get("TAXON_IDX_MAPPING_FILE_PATH")
        results_base_dir = os.environ.get("RESULTS_BASE_DIR")

        # Without env variables being absolute paths, application won't start
        if cluster_f is None or not os.path.isabs(cluster_f):
            logger.error("[ERROR] CLUSTER_FILE_PATH should be an absolute path.")
            sys.exit(1)
        if sequence_ids_f is None or not os.path.isabs(sequence_ids_f):
            logger.error("[ERROR] SEQUENCE_IDS_FILE_PATH should be an absolute path.")
            sys.exit(1)
        if taxon_idx_mapping_file is None or not os.path.isabs(taxon_idx_mapping_file):
            logger.error("[ERROR] TAXON_IDX_MAPPING_FILE_PATH should be an absolute path.")  # fmt:skip
            sys.exit(1)
        if results_base_dir is None or not os.path.isabs(results_base_dir):
            logger.error("[ERROR] RESULTS_BASE_DIR should be an absolute path.")
            sys.exit(1)

        try:
            check_file(cluster_f, install_kinfin=True)
            check_file(sequence_ids_f, install_kinfin=True)
            check_file(taxon_idx_mapping_file, install_kinfin=True)
        except FileNotFoundError as e:
            logger.error(e)
            sys.exit(1)

        run_server(
            args=args,
            nodesdb_f=nodesdb_f,
            go_mapping_f=go_mapping_f,
            ipr_mapping_f=ipr_mapping_f,
            pfam_mapping_f=pfam_mapping_f,
            cluster_f=cluster_f,
            sequence_ids_f=sequence_ids_f,
            results_base_dir=results_base_dir,
            taxon_idx_mapping_file=taxon_idx_mapping_file,
        )
    elif isinstance(args, InputData):
        # run the cli script
        run_cli(args)

    else:
        logger.error("[ERROR] - invalid input provided.")
        sys.exit(1)
