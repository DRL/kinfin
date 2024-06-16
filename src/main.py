import os
import sys

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
        pass
    elif isinstance(args, InputData):
        # run the cli script
        pass

    else:
        print("[ERROR] - invalid input provided.")
        sys.exit(1)
