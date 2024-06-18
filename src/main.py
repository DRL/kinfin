import os
import sys

from api import run_server
from cli import run_cli
from cli.commands import parse_args
from core.input import InputData, ServeArgs
from core.utils import check_file

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
        print(e)
        sys.exit(1)

    args = parse_args(nodesdb_f, pfam_mapping_f, ipr_mapping_f, go_mapping_f)

    if isinstance(args, ServeArgs):
        # run the api server
        cluster_f = os.path.join(base_dir, "data/api/advanced/Orthogroups.txt")
        sequence_ids_f = os.path.join(base_dir, "data/api/advanced/SequenceIDs.txt")
        taxon_idx_mapping_file = os.path.join(base_dir, "data/api/advanced/taxon_idx_mapping.json")  # fmt:skip
        try:
            check_file(cluster_f, install_kinfin=True)
            check_file(sequence_ids_f, install_kinfin=True)
            check_file(taxon_idx_mapping_file, install_kinfin=True)
        except FileNotFoundError as e:
            print(e)
            sys.exit(1)

        run_server(
            args=args,
            base_dir=base_dir,
            cluster_f=cluster_f,
            go_mapping_f=go_mapping_f,
            ipr_mapping_f=ipr_mapping_f,
            nodesdb_f=nodesdb_f,
            pfam_mapping_f=pfam_mapping_f,
            sequence_ids_f=sequence_ids_f,
            taxon_idx_mapping_file=taxon_idx_mapping_file,
        )
    elif isinstance(args, InputData):
        # run the cli script
        run_cli(args)

    else:
        print("[ERROR] - invalid input provided.")
        sys.exit(1)
