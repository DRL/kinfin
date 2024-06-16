import os
import shutil
import time
from typing import Dict, List, Literal, Optional

from core.alo_collections import AloCollection
from core.build import (
    build_AloCollection,
    build_AloCollection_from_json,
    build_ClusterCollection,
    build_ProteinCollection,
)
from core.clusters import Cluster, ClusterCollection
from core.input import InputData
from core.logic import get_ALO_cluster_cardinality, get_attribute_cluster_type
from core.proteins import ProteinCollection
from core.utils import median, progress, statistic


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

    def analyse_cluster(self, cluster: Cluster):
        """
        Analyzes the given cluster and updates various attributes and statistics related to
        different levels and attributes in the ALOCollection.

        Args:
            cluster (Cluster): The cluster to analyze.

        Notes:
            This method performs analysis on the provided cluster against the ALOCollection's
            tree structure. It computes various statistics such as protein counts, proteome
            coverage, and determines cluster types for each attribute level.

            The analysis results are stored in the `cluster` object itself and include:
            - `protein_counts_of_proteomes_by_level_by_attribute`
            - `proteome_coverage_by_level_by_attribute`
            - `implicit_protein_ids_by_proteome_id_by_level_by_attribute`
            - `cluster_type_by_attribute`
            - `protein_median`

            Additionally, the method updates the ALOCollection's ALO objects with cluster
            information including protein IDs, statistics, and test results where applicable.

        Returns:
            None

        Raises:
            Any exceptions that may occur during the execution of statistical tests or data
            retrieval operations.

        """
        protein_get_by_proteome_id = cluster.protein_ids_by_proteome_id.get

        implicit_protein_ids_by_proteome_id_by_level_by_attribute: Dict[str, Dict[str, Dict[str, List[str]]]] = {}  # fmt:skip
        cluster_type_by_attribute: Dict[str, Literal["singleton", "shared", "specific"]] = {}  # fmt:skip
        protein_counts_of_proteomes_by_level_by_attribute:Dict[str,Dict[str,List[int]]] = {}  # fmt:skip
        proteome_coverage_by_level_by_attribute: Dict[str,Dict[str,float]] = {}  # fmt:skip

        if self.aloCollection.tree_ete:
            for node in self.aloCollection.tree_ete.traverse("levelorder"):  # type: ignore
                intersection = cluster.proteome_ids.intersection(node.proteome_ids)  # type: ignore
                difference = cluster.proteome_ids.difference(node.proteome_ids)  # type: ignore
                if len(intersection) == 0:
                    # Nothing to see here ...
                    node.counts["absent"] += 1  # type: ignore
                else:
                    if cluster.singleton is True:
                        # This is a singleton
                        node.counts["singleton"] += 1  # type: ignore
                        node.apomorphic_cluster_counts["singletons"] += 1  # type: ignore
                    elif len(difference) > 0:
                        # This is a 'shared' cluster
                        node.counts["shared"] += 1  # type: ignore
                    elif len(difference) == 0:
                        # This is a node 'specific' cluster
                        node.counts["specific"] += 1  # type: ignore
                        if cluster.proteome_count == 1:
                            # But it only belongs to one proteome
                            node.apomorphic_cluster_counts["non_singletons"] += 1  # type: ignore
                        else:
                            # It has more than one proteome
                            child_nodes_covered = []
                            child_node_proteome_coverage_strings = []
                            child_node_proteome_ids_covered_count = 0
                            for child_node in node.get_children():
                                if child_node.proteome_ids.isdisjoint(
                                    cluster.proteome_ids
                                ):
                                    # No child node proteomes are not in cluster
                                    child_nodes_covered.append(False)
                                else:
                                    # At least on child node proteome in cluster
                                    child_nodes_covered.append(True)
                                    child_node_proteome_ids_covered_count = len(
                                        cluster.proteome_ids.intersection(
                                            child_node.proteome_ids
                                        )
                                    )
                                    child_node_proteome_coverage_strings.append(
                                        f"{child_node.name}=({child_node_proteome_ids_covered_count}/{len(child_node.proteome_ids)})"
                                    )
                            if all(child_nodes_covered):
                                # At least one proteome of each child node in cluster
                                # => SYNAPOMORPHY
                                node_proteome_coverage = len(intersection) / len(
                                    node.proteome_ids  # type: ignore
                                )
                                node_cluster_type = ""
                                if node_proteome_coverage == 1.0:
                                    node_cluster_type = "complete_presence"
                                else:
                                    node_cluster_type = "partial_absence"
                                node.synapomorphic_cluster_counts[node_cluster_type] += 1  # type: ignore

                                node.synapomorphic_cluster_strings.append(  # type: ignore
                                    (
                                        cluster.cluster_id,
                                        node.name,
                                        node_cluster_type,
                                        "{0:.3}".format(node_proteome_coverage),
                                        ";".join(child_node_proteome_coverage_strings),
                                        ",".join(sorted(intersection)),
                                    )
                                )
        for attribute in self.aloCollection.attributes:
            protein_counts_of_proteomes_by_level_by_attribute[attribute] = {}
            proteome_coverage_by_level_by_attribute[attribute] = {}
            implicit_protein_ids_by_proteome_id_by_level_by_attribute[attribute] = {}
            protein_ids_by_level: Dict[str, List[str]] = {}
            protein_length_stats_by_level: Dict[str, Dict[str, int | float]] = {}
            explicit_protein_count_by_proteome_id_by_level: Dict[str, Dict[str,int]] = {}  # fmt:skip

            for level in self.aloCollection.ALO_by_level_by_attribute[attribute]:
                protein_ids_by_proteome_id: Dict[str, List[str]] = {}
                protein_count_by_proteome_id: Dict[str, int] = {}
                protein_ids_by_level[level] = []

                ALO = self.aloCollection.ALO_by_level_by_attribute[attribute][level]

                if ALO is not None:
                    for proteome_id in ALO.proteomes_list:
                        protein_ids: List[str] = list(protein_get_by_proteome_id(proteome_id, []))  # fmt:skip
                        protein_ids_by_level[level] += protein_ids
                        protein_count_by_proteome_id[proteome_id] = len(protein_ids)
                        if not protein_count_by_proteome_id[proteome_id] == 0:
                            protein_ids_by_proteome_id[proteome_id] = protein_ids

                    if protein_ids_by_proteome_id:
                        implicit_protein_ids_by_proteome_id_by_level_by_attribute[attribute][level] = protein_ids_by_proteome_id  # fmt:skip

                    explicit_protein_count_by_proteome_id_by_level[level] = protein_count_by_proteome_id  # fmt:skip

                    protein_length_stats_by_level[level] = self.proteinCollection.get_protein_length_stats(protein_ids_by_level[level])  # fmt:skip

                    protein_counts_of_proteomes_by_level_by_attribute[attribute][level] = [
                        protein_count
                        for _, protein_count in list(
                            protein_count_by_proteome_id.items()
                        )
                    ]  # fmt:skip

                cluster_type_by_attribute[attribute] = get_attribute_cluster_type(
                    cluster.singleton,
                    implicit_protein_ids_by_proteome_id_by_level_by_attribute[attribute]  # fmt:skip
                )

            for level in self.aloCollection.ALO_by_level_by_attribute[attribute]:
                ALO = self.aloCollection.ALO_by_level_by_attribute[attribute][level]
                if ALO is not None:
                    proteome_coverage_by_level_by_attribute[attribute][level] = (
                        len(
                            implicit_protein_ids_by_proteome_id_by_level_by_attribute[
                                attribute
                            ].get(level, [])
                        )
                        / ALO.proteome_count
                    )

                    ALO_cluster_status: str = ""
                    ALO_cluster_cardinality: Optional[str] = None
                    mwu_pvalue: Optional[float] = None
                    mwu_log2_mean: Optional[float] = None
                    mean_ALO_count: Optional[float] = None
                    mean_non_ALO_count: Optional[float] = None

                    if level not in implicit_protein_ids_by_proteome_id_by_level_by_attribute[attribute]:  # fmt:skip
                        ALO_cluster_status = "absent"
                    else:
                        ALO_cluster_status = "present"

                        if not cluster_type_by_attribute[attribute] == "singleton":
                            ALO_proteome_counts_in_cluster: List[int] = [
                                count
                                for _, count in list(
                                    explicit_protein_count_by_proteome_id_by_level[
                                        level
                                    ].items()
                                )
                            ]
                            ALO_cluster_cardinality = get_ALO_cluster_cardinality(
                                ALO_proteome_counts_in_cluster=ALO_proteome_counts_in_cluster,
                                fuzzy_count=self.inputData.fuzzy_count,
                                fuzzy_fraction=self.inputData.fuzzy_fraction,
                                fuzzy_range=self.inputData.fuzzy_range,
                            )
                            if cluster_type_by_attribute[attribute] == "shared":
                                non_ALO_levels = [
                                    non_ALO_level
                                    for non_ALO_level in explicit_protein_count_by_proteome_id_by_level
                                    if not non_ALO_level == level
                                ]
                                non_ALO_proteome_counts_in_cluster = []
                                for non_ALO_level in non_ALO_levels:
                                    for (
                                        proteome_id
                                    ) in explicit_protein_count_by_proteome_id_by_level[
                                        non_ALO_level
                                    ]:
                                        non_ALO_proteome_counts_in_cluster.append(
                                            explicit_protein_count_by_proteome_id_by_level[
                                                non_ALO_level
                                            ][
                                                proteome_id
                                            ]
                                        )
                                (
                                    mwu_pvalue,
                                    mwu_log2_mean,
                                    mean_ALO_count,
                                    mean_non_ALO_count,
                                ) = statistic(
                                    count_1=ALO_proteome_counts_in_cluster,
                                    count_2=non_ALO_proteome_counts_in_cluster,
                                    test=self.inputData.test,
                                    min_proteomes=self.inputData.min_proteomes,
                                )

                    ALO.add_cluster(
                        cluster=cluster,
                        attribute_cluster_type=cluster_type_by_attribute[attribute],
                        ALO_cluster_status=ALO_cluster_status,
                        ALO_protein_length_stats=protein_length_stats_by_level[level],
                        ALO_protein_ids_in_cluster=protein_ids_by_level[level],
                        ALO_cluster_cardinality=ALO_cluster_cardinality,
                        mwu_pvalue=mwu_pvalue,
                        mwu_log2_mean=mwu_log2_mean,
                        mean_ALO_count=mean_ALO_count,
                        mean_non_ALO_count=mean_non_ALO_count,
                    )

        cluster.protein_counts_of_proteomes_by_level_by_attribute = protein_counts_of_proteomes_by_level_by_attribute  # fmt:skip
        cluster.proteome_coverage_by_level_by_attribute = proteome_coverage_by_level_by_attribute  # fmt:skip
        cluster.implicit_protein_ids_by_proteome_id_by_level_by_attribute = implicit_protein_ids_by_proteome_id_by_level_by_attribute  # fmt:skip
        cluster.cluster_type_by_attribute = cluster_type_by_attribute  # fmt:skip
        cluster.protein_median = median(
            [
                count
                for count in protein_counts_of_proteomes_by_level_by_attribute["all"][
                    "all"
                ]
                if not count == 0
            ]
        )

    def analyse_clusters(self) -> None:
        if self.clusterCollection.inferred_singletons_count:
            print(f"[STATUS]\t - Clusters found = {self.clusterCollection.cluster_count} (of which {self.clusterCollection.inferred_singletons_count} were inferred singletons)")  # fmt:skip

        else:
            print(f"[STATUS]\t - Clusters found = {self.clusterCollection.cluster_count}")  # fmt:skip

        parse_steps = self.clusterCollection.cluster_count / 100

        print("[STATUS] - Analysing clusters ...")
        analyse_clusters_start = time.time()
        for idx, cluster in enumerate(self.clusterCollection.cluster_list):
            self.analyse_cluster(cluster)
            progress(idx + 1, parse_steps, self.clusterCollection.cluster_count)
        analyse_clusters_end = time.time()
        analyse_clusters_elapsed = analyse_clusters_end - analyse_clusters_start
        print(f"[STATUS] - Took {analyse_clusters_elapsed}s to analyse clusters")
