import os
import shutil
import time
from collections import Counter
from typing import Any, Dict, Generator, List, Literal, Optional

import matplotlib as mat
import numpy as np

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

mat.use("agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

plt.style.use("ggplot")
mat.rc("ytick", labelsize=20)
mat.rc("xtick", labelsize=20)
axis_font = {"size": "20"}
mat.rcParams.update({"font.size": 22})


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

    def plot_cluster_sizes(self):
        """
        Plot the distribution of cluster sizes based on the protein counts in each cluster.

        This method generates a scatter plot showing the distribution of cluster sizes
        based on the number of proteins in each cluster.

        Args:
        - None

        Returns:
        - None
        """
        cluster_protein_count = []
        for cluster in self.clusterCollection.cluster_list:
            cluster_protein_count.append(cluster.protein_count)
        cluster_protein_counter = Counter(cluster_protein_count)
        count_plot_f = os.path.join(
            self.dirs["main"], f"cluster_size_distribution.{self.inputData.plot_format}"
        )
        f, ax = plt.subplots(figsize=self.inputData.plotsize)
        ax.set_facecolor("white")
        x_values = []
        y_values = []
        for value, count in list(cluster_protein_counter.items()):
            x_values.append(value)
            y_values.append(count)
        x_array = np.array(x_values)
        y_array = np.array(y_values)
        ax.scatter(x_array, y_array, marker="o", alpha=0.8, s=100)
        ax.set_xlabel("Cluster size", fontsize=self.inputData.fontsize)
        ax.set_ylabel("Count", fontsize=self.inputData.fontsize)
        ax.set_yscale("log")
        ax.set_xscale("log")
        plt.margins(0.8)
        plt.gca().set_ylim(bottom=0.8)
        plt.gca().set_xlim(left=0.8)
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
        f.tight_layout()

        ax.grid(True, linewidth=1, which="major", color="lightgrey")
        ax.grid(True, linewidth=0.5, which="minor", color="lightgrey")
        print(f"[STATUS] - Plotting {count_plot_f}")
        f.savefig(count_plot_f, format=self.inputData.plot_format)
        plt.close()

    def get_header_line(self, filetype, attribute):
        """
        Generate a header line based on the specified filetype and attribute.

        Args:
        - filetype (str): The type of file for which the header line is generated.
                          Possible values include "attribute_metrics", "cluster_metrics_ALO",
                          "cluster_metrics", "cluster_metrics_domains", "cluster_metrics_domains_detailed",
                          "cafe", "pairwise_representation_test", "cluster_1to1s_ALO".
        - attribute (str): The specific attribute or level associated with the header line.

        Returns:
        - str: The generated header line as a tab-separated string.

        Raises:
        - ValueError: If the provided filetype is not recognized.

        """
        if filetype == "attribute_metrics":
            attribute_metrics_header = []
            attribute_metrics_header.append("#attribute")
            attribute_metrics_header.append("taxon_set")
            attribute_metrics_header.append("cluster_total_count")
            attribute_metrics_header.append("protein_total_count")
            attribute_metrics_header.append("protein_total_span")
            attribute_metrics_header.append("singleton_cluster_count")
            attribute_metrics_header.append("singleton_protein_count")
            attribute_metrics_header.append("singleton_protein_span")
            attribute_metrics_header.append("specific_cluster_count")
            attribute_metrics_header.append("specific_protein_count")
            attribute_metrics_header.append("specific_protein_span")
            attribute_metrics_header.append("shared_cluster_count")
            attribute_metrics_header.append("shared_protein_count")
            attribute_metrics_header.append("shared_protein_span")
            attribute_metrics_header.append("specific_cluster_true_1to1_count")
            attribute_metrics_header.append("specific_cluster_fuzzy_count")
            attribute_metrics_header.append("shared_cluster_true_1to1_count")
            attribute_metrics_header.append("shared_cluster_fuzzy_count")
            attribute_metrics_header.append("absent_cluster_total_count")
            attribute_metrics_header.append("absent_cluster_singleton_count")
            attribute_metrics_header.append("absent_cluster_specific_count")
            attribute_metrics_header.append("absent_cluster_shared_count")
            attribute_metrics_header.append("TAXON_count")
            attribute_metrics_header.append("TAXON_taxa")
            return "\t".join(attribute_metrics_header)
        elif filetype == "cluster_metrics_ALO":
            cluster_metrics_ALO_header = []
            cluster_metrics_ALO_header.append("#cluster_id")
            cluster_metrics_ALO_header.append("cluster_status")
            cluster_metrics_ALO_header.append("cluster_type")
            cluster_metrics_ALO_header.append("cluster_protein_count")
            cluster_metrics_ALO_header.append("cluster_proteome_count")
            cluster_metrics_ALO_header.append("TAXON_protein_count")
            cluster_metrics_ALO_header.append("TAXON_mean_count")
            cluster_metrics_ALO_header.append("non_taxon_mean_count")
            cluster_metrics_ALO_header.append("representation")
            cluster_metrics_ALO_header.append("log2_mean(TAXON/others)")
            cluster_metrics_ALO_header.append("pvalue(TAXON vs. others)")
            cluster_metrics_ALO_header.append("TAXON_coverage")
            cluster_metrics_ALO_header.append("TAXON_count")
            cluster_metrics_ALO_header.append("non_TAXON_count")
            cluster_metrics_ALO_header.append("TAXON_taxa")
            cluster_metrics_ALO_header.append("non_TAXON_taxa")
            # for domain_source in clusterCollection.domain_sources:
            #    cluster_metrics_ALO_header.append(domain_source)
            return "\t".join(cluster_metrics_ALO_header)
        elif filetype == "cluster_metrics":
            cluster_metrics_header = []
            cluster_metrics_header.append("#cluster_id")
            cluster_metrics_header.append("cluster_protein_count")
            cluster_metrics_header.append("protein_median_count")
            cluster_metrics_header.append("TAXON_count")
            cluster_metrics_header.append("attribute")
            cluster_metrics_header.append("attribute_cluster_type")
            cluster_metrics_header.append("protein_span_mean")
            cluster_metrics_header.append("protein_span_sd")
            cluster_metrics_header += [
                "%s_count" % level
                for level in sorted(
                    self.aloCollection.ALO_by_level_by_attribute[attribute]
                )
            ]
            if not attribute == "TAXON":
                cluster_metrics_header += [
                    "%s_median" % level
                    for level in sorted(
                        self.aloCollection.ALO_by_level_by_attribute[attribute]
                    )
                ]
                cluster_metrics_header += [
                    "%s_cov" % level
                    for level in sorted(
                        self.aloCollection.ALO_by_level_by_attribute[attribute]
                    )
                ]
            return "\t".join(cluster_metrics_header)
        elif filetype == "cluster_metrics_domains":
            cluster_metrics_domains_header = []
            cluster_metrics_domains_header.append("#cluster_id")
            cluster_metrics_domains_header.append("cluster_protein_count")
            cluster_metrics_domains_header.append("TAXON_count")
            cluster_metrics_domains_header.append("protein_span_mean")
            cluster_metrics_domains_header.append("protein_span_sd")
            cluster_metrics_domains_header.append("fraction_secreted")
            for domain_source in self.clusterCollection.domain_sources:
                cluster_metrics_domains_header.append(domain_source)
                cluster_metrics_domains_header.append(f"{domain_source}_entropy")
            return "\t".join(cluster_metrics_domains_header)
        elif filetype == "cluster_metrics_domains_detailed":
            cluster_metrics_domains_detailed_header = []
            cluster_metrics_domains_detailed_header.append("#cluster_id")
            cluster_metrics_domains_detailed_header.append("domain_source")
            cluster_metrics_domains_detailed_header.append("domain_id")
            cluster_metrics_domains_detailed_header.append("domain_description")
            cluster_metrics_domains_detailed_header.append("protein_count")
            cluster_metrics_domains_detailed_header.append("protein_count_with_domain")
            cluster_metrics_domains_detailed_header.append("TAXA_with_domain_fraction")
            cluster_metrics_domains_detailed_header.append("TAXA_with_domain")
            cluster_metrics_domains_detailed_header.append("TAXA_without_domain")
            return "\t".join(cluster_metrics_domains_detailed_header)
        elif filetype == "cafe":
            cafe_header = []
            cafe_header.append("#ID")
            for level in sorted(self.aloCollection.ALO_by_level_by_attribute["TAXON"]):
                cafe_header.append(level)
            return "\t".join(cafe_header)
        elif filetype == "pairwise_representation_test":
            pairwise_representation_test_header = []
            pairwise_representation_test_header.append("#cluster_id")
            pairwise_representation_test_header.append("TAXON_1")
            pairwise_representation_test_header.append("TAXON_1_mean")
            pairwise_representation_test_header.append("TAXON_2")
            pairwise_representation_test_header.append("TAXON_2_mean")
            pairwise_representation_test_header.append("log2_mean(TAXON_1/TAXON_2)")
            pairwise_representation_test_header.append(
                "mwu_pvalue(TAXON_1 vs. TAXON_2)"
            )
            # pairwise_representation_test_header.append("go_terms")
            # for domain_source in clusterCollection.domain_sources:
            #    pairwise_representation_test_header.append(domain_source)
            return "\t".join(pairwise_representation_test_header)
        elif filetype == "cluster_1to1s_ALO":
            cluster_1to1s_ALO_header = []
            cluster_1to1s_ALO_header.append("#cluster_id")
            cluster_1to1s_ALO_header.append("cluster_type")
            cluster_1to1s_ALO_header.append("1to1_type")
            cluster_1to1s_ALO_header.append("proteome_count")
            cluster_1to1s_ALO_header.append("percentage_at_target_count")
            return "\t".join(cluster_1to1s_ALO_header)
        else:
            error_msg = f"[ERROR] {filetype} is not a valid header 'filetype'"
            raise ValueError(error_msg)

    def pairwise_representation_test(
        self, clusterObj, attribute, level, levels_seen, levels
    ) -> Generator[Any, Any, Any]:
        """
        Conducts pairwise statistical tests between the protein counts of a given cluster
        at a specific level (`level`) and all other levels (`other_level`) within a set of `levels`.

        Parameters:
        - clusterObj: The cluster object for which pairwise tests are conducted.
        - attribute: The attribute or level associated with the cluster metrics.
        - level: The specific level within the attribute for comparisons.
        - levels_seen: A set containing levels that have already been processed.
        - levels: A collection of all levels within the attribute for potential comparisons.

        Yields:
        - List: A list containing the following statistical results:
        [cluster_id, level, other_level, mean_ALO_count, mean_non_ALO_count, mwu_log2_m

        """

        for other_level in set(levels).difference(levels_seen):
            if not other_level == level:
                other_ALO = self.aloCollection.ALO_by_level_by_attribute[attribute][
                    other_level
                ]
                if (
                    other_ALO
                    and len(clusterObj.proteome_ids.intersection(other_ALO.proteomes))
                    >= 2
                ):
                    protein_counts_level = [
                        count
                        for count in clusterObj.protein_counts_of_proteomes_by_level_by_attribute[
                            attribute
                        ][
                            level
                        ]
                        if count > 0
                    ]
                    protein_counts_other_level = [
                        count
                        for count in clusterObj.protein_counts_of_proteomes_by_level_by_attribute[
                            attribute
                        ][
                            other_level
                        ]
                        if count > 0
                    ]
                    if protein_counts_level and protein_counts_other_level:
                        (
                            mwu_pvalue,
                            mwu_log2_mean,
                            mean_ALO_count,
                            mean_non_ALO_count,
                        ) = statistic(
                            protein_counts_level,
                            protein_counts_other_level,
                            self.inputData.test,
                            self.inputData.min_proteomes,
                        )
                        yield [
                            clusterObj.cluster_id,
                            level,
                            other_level,
                            mean_ALO_count,
                            mean_non_ALO_count,
                            mwu_log2_mean,
                            mwu_pvalue,
                        ]
                        # pvalue = None
                        # try:
                        #     pvalue = scipy.stats.mannwhitneyu(protein_counts_level, protein_counts_other_level, alternative="two-sided")[1]
                        # except:
                        #     pvalue = 1.0
                        # mean_level = mean(protein_counts_level)
                        # mean_other_level = mean(protein_counts_other_level)
                        # log2fc_mean = log((mean_level/mean_other_level), 2)
                        # yield [clusterObj.cluster_id, level, other_level, mean_level, mean_other_level, log2fc_mean, pvalue]
