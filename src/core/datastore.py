import os
import shutil
from typing import Optional

from core.alo_collections import AloCollection
from core.build import (
    build_AloCollection,
    build_AloCollection_from_json,
    build_ClusterCollection,
    build_ProteinCollection,
)
from core.clusters import ClusterCollection
from core.input import InputData
from core.proteins import ProteinCollection


class DataFactory:
    def __init__(self, inputData: InputData) -> None:
        self.dirs = {}
        self.inputData: InputData = inputData
        if isinstance(self.inputData.config_data, str):
            self.aloCollection: AloCollection = build_AloCollection(
                config_f=self.inputData.config_data,
                nodesdb_f=self.inputData.nodesdb_f,
                tree_f=self.inputData.tree_f,
                taxranks=self.inputData.taxranks,
            )
        elif self.inputData.taxon_idx_mapping_file is not None:
            self.aloCollection: AloCollection = build_AloCollection_from_json(
                nodesdb_f=self.inputData.nodesdb_f,
                tree_f=self.inputData.tree_f,
                taxranks=self.inputData.taxranks,
                json_list=self.inputData.config_data,
                taxon_idx_mapping_file=self.inputData.taxon_idx_mapping_file,
            )
            pass
        else:
            raise ValueError("[ERROR] - Either provide config file or json")

        self.proteinCollection: ProteinCollection = build_ProteinCollection(
            aloCollection=self.aloCollection,
            fasta_dir=self.inputData.fasta_dir,
            go_mapping_f=self.inputData.go_mapping_f,
            functional_annotation_f=self.inputData.functional_annotation_f,
            ipr_mapping=self.inputData.ipr_mapping,
            ipr_mapping_f=self.inputData.ipr_mapping_f,
            pfam_mapping=self.inputData.pfam_mapping,
            pfam_mapping_f=self.inputData.pfam_mapping_f,
            sequence_ids_f=self.inputData.sequence_ids_f,
            species_ids_f=self.inputData.species_ids_f,
        )
        self.clusterCollection: ClusterCollection = build_ClusterCollection(
            cluster_f=self.inputData.cluster_f,
            proteinCollection=self.proteinCollection,
            infer_singletons=self.inputData.infer_singletons,
        )

    def setup_dirs(self):
        """
        Set up output directories for storing results and attributes.
        """
        output_path: Optional[str] = self.inputData.output_path

        if output_path:
            if not os.path.isabs(output_path):
                output_path = os.path.abspath(output_path)
        else:
            output_path = os.path.join(os.getcwd(), "kinfin_results")

        self.dirs["main"] = output_path
        print(f"[STATUS] - Output directories in \n\t{output_path}")
        if os.path.exists(output_path):
            print("[STATUS] - Directory exists. Deleting directory ...")
            shutil.rmtree(output_path)

        print("[STATUS] - Creating directories ...")
        os.mkdir(output_path)
        for attribute in self.aloCollection.attributes:
            print(attribute)
            attribute_path = os.path.join(output_path, attribute)
            self.dirs[attribute] = attribute_path
            if not os.path.exists(attribute_path):
                print(f"\t{attribute_path}")
                os.mkdir(attribute_path)

        if self.aloCollection.tree_ete is not None:
            tree_path = os.path.join(output_path, "tree")
            node_chart_path = os.path.join(tree_path, "charts")
            node_header_path = os.path.join(tree_path, "headers")

            if not os.path.exists(tree_path):
                print(f"\t{tree_path}")
                os.mkdir(tree_path)
                self.dirs["tree"] = tree_path

                print(f"\t{node_chart_path}")
                os.mkdir(node_chart_path)
                self.dirs["tree_charts"] = node_chart_path

                if self.inputData.plot_tree:
                    print(f"\t{node_header_path}")
                    os.mkdir(node_header_path)
                    self.dirs["tree_headers"] = node_header_path
