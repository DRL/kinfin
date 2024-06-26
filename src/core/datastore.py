import os
import shutil
import time
from collections import Counter
from typing import Any, Dict, Generator, List, Literal, Optional, Set, Tuple

import matplotlib as mat
import numpy as np
from matplotlib.lines import Line2D

from core.alo import AttributeLevel
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
from core.utils import logger, median, progress, statistic

mat.use("agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, NullFormatter

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
        logger.info(f"[STATUS] - Output directories in")
        logger.info(f"\t{output_path}")
        if os.path.exists(output_path):
            logger.info("[STATUS] - Directory exists. Deleting directory ...")
            shutil.rmtree(output_path)

        logger.info("[STATUS] - Creating directories ...")
        os.mkdir(output_path)
        for attribute in self.aloCollection.attributes:
            attribute_path = os.path.join(output_path, attribute)
            self.dirs[attribute] = attribute_path
            if not os.path.exists(attribute_path):
                logger.info(f"\t{attribute_path}")
                os.mkdir(attribute_path)

        if self.aloCollection.tree_ete is not None:
            tree_path = os.path.join(output_path, "tree")
            node_chart_path = os.path.join(tree_path, "charts")
            node_header_path = os.path.join(tree_path, "headers")

            if not os.path.exists(tree_path):
                logger.info(f"\t{tree_path}")
                os.mkdir(tree_path)
                self.dirs["tree"] = tree_path

                logger.info(f"\t{node_chart_path}")
                os.mkdir(node_chart_path)
                self.dirs["tree_charts"] = node_chart_path

                if self.inputData.plot_tree:
                    logger.info(f"\t{node_header_path}")
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
            logger.info(f"[STATUS]\t - Clusters found = {self.clusterCollection.cluster_count} (of which {self.clusterCollection.inferred_singletons_count} were inferred singletons)")  # fmt:skip

        else:
            logger.info(f"[STATUS]\t - Clusters found = {self.clusterCollection.cluster_count}")  # fmt:skip

        parse_steps = self.clusterCollection.cluster_count / 100

        logger.info("[STATUS] - Analysing clusters ...")
        analyse_clusters_start = time.time()
        for idx, cluster in enumerate(self.clusterCollection.cluster_list):
            self.analyse_cluster(cluster)
            progress(idx + 1, parse_steps, self.clusterCollection.cluster_count)
        analyse_clusters_end = time.time()
        analyse_clusters_elapsed = analyse_clusters_end - analyse_clusters_start
        logger.info(f"[STATUS] - Took {analyse_clusters_elapsed}s to analyse clusters")

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
        logger.info(f"[STATUS] - Plotting {count_plot_f}")
        f.savefig(count_plot_f, format=self.inputData.plot_format)
        plt.close()

    def get_header_line(self, filetype: str, attribute: str) -> str:
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
        self,
        cluster: Cluster,
        attribute: str,
        level: str,
        levels_seen: Set[str],
        levels: List[str],
    ) -> Generator[Any, Any, Any]:
        """
        Conducts pairwise statistical tests between the protein counts of a given cluster
        at a specific level (`level`) and all other levels (`other_level`) within a set of `levels`.

        Parameters:
        - cluster: The cluster object for which pairwise tests are conducted.
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
                    and len(cluster.proteome_ids.intersection(other_ALO.proteomes)) >= 2
                ):
                    protein_counts_level = [
                        count
                        for count in cluster.protein_counts_of_proteomes_by_level_by_attribute[
                            attribute
                        ][
                            level
                        ]
                        if count > 0
                    ]
                    protein_counts_other_level = [
                        count
                        for count in cluster.protein_counts_of_proteomes_by_level_by_attribute[
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
                            cluster.cluster_id,
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
                        # yield [cluster.cluster_id, level, other_level, mean_level, mean_other_level, log2fc_mean, pvalue]

    def get_attribute_metrics(self, ALO: AttributeLevel) -> str:
        """
        Retrieves various metrics related to an Attribute-Level Orthology (ALO) object.

        Parameters:
        - ALO: An Attribute-Level Orthology (ALO) object containing cluster and protein metrics.

        Returns:
        - str: A tab-separated string containing the following metrics:
        [attribute, level, cluster_total_count, protein_total_count, protein_total_span,
        singleton_cluster_count, singleton_protein_count, singleton_protein_span,
        specific_cluster_count, specific_protein_count, specific_protein_span,
        shared_cluster_count, shared_protein_count, shared_protein_span,
        specific_cluster_true_1to1_count, specific_cluster_fuzzy_count,
        shared_cluster_true_1to1_count, shared_cluster_fuzzy_count,
        absent_cluster_total_count, absent_cluster_singleton_count,
        absent_cluster_specific_count, absent_cluster_shared_count,
        proteome_count, TAXON_taxa]
        """
        attribute_metrics = []
        attribute_metrics.append(ALO.attribute)
        attribute_metrics.append(ALO.level)
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_status_by_cluster_type("present", "total")
        )
        attribute_metrics.append(ALO.get_protein_count_by_cluster_type("total"))
        attribute_metrics.append(ALO.get_protein_span_by_cluster_type("total"))
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_status_by_cluster_type(
                "present", "singleton"
            )
        )
        attribute_metrics.append(ALO.get_protein_count_by_cluster_type("singleton"))
        attribute_metrics.append(ALO.get_protein_span_by_cluster_type("singleton"))
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_status_by_cluster_type(
                "present", "specific"
            )
        )
        attribute_metrics.append(ALO.get_protein_count_by_cluster_type("specific"))
        attribute_metrics.append(ALO.get_protein_span_by_cluster_type("specific"))
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_status_by_cluster_type("present", "shared")
        )
        attribute_metrics.append(ALO.get_protein_count_by_cluster_type("shared"))
        attribute_metrics.append(ALO.get_protein_span_by_cluster_type("shared"))
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_cardinality_by_cluster_type(
                "specific", "true"
            )
        )
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_cardinality_by_cluster_type(
                "specific", "fuzzy"
            )
        )
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_cardinality_by_cluster_type(
                "shared", "true"
            )
        )
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_cardinality_by_cluster_type(
                "shared", "fuzzy"
            )
        )
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_status_by_cluster_type("absent", "total")
        )
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_status_by_cluster_type(
                "absent", "singleton"
            )
        )
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_status_by_cluster_type(
                "absent", "specific"
            )
        )
        attribute_metrics.append(
            ALO.get_cluster_count_by_cluster_status_by_cluster_type("absent", "shared")
        )
        attribute_metrics.append(ALO.proteome_count)
        attribute_metrics.append(ALO.get_proteomes())
        return "\t".join([str(field) for field in attribute_metrics])

    def plot_count_comparisons_volcano(
        self,
        pairwise_representation_test_by_pair_by_attribute,
    ) -> None:
        """
        Plot volcano plots for pairwise comparisons based on statistical tests.

        Parameters:
        - pairwise_representation_test_by_pair_by_attribute (dict): A dictionary where keys are attributes
          and values are dictionaries mapping pairs to their test results.

        Notes:
        - This function generates volcano plots for each pair of attributes, showing the relationship between
          log2 fold change and p-values.

        Returns:
        - None: Plots are saved as files.
        """
        # [cluster.cluster_id, level, other_level, mean_level, mean_other_level, log2fc_mean, pvalue]
        for attribute in pairwise_representation_test_by_pair_by_attribute:
            for pair in pairwise_representation_test_by_pair_by_attribute[attribute]:
                pair_list = list(pair)
                x_label = pair_list[0]
                y_label = pair_list[1]
                pair_data = pairwise_representation_test_by_pair_by_attribute[
                    attribute
                ][pair]
                pair_data_count = len(pair_data)
                p_values = []
                log2fc_values = []
                for data in pair_data:
                    log2fc_values.append(data[5])
                    pvalue = data[6]
                    if pvalue == 0.0:
                        pvalue = 0.01 / (
                            pair_data_count + 1
                        )  # if value is 0.0 it gets set beyond the bonferroni corrected 0.01
                    p_values.append(pvalue)
                if p_values:
                    pairwise_representation_test_f = os.path.join(
                        self.dirs[attribute],
                        f"{attribute}.pairwise_representation_test.{"_".join(pair_list)}.{self.inputData.plot_format}",
                    )

                    plt.figure(1, figsize=self.inputData.plotsize)
                    # f, ax = plt.subplots(figsize=self.inputData.plotsize)
                    # ax.set_facecolor('white')
                    # definitions for the axes
                    left, width = 0.1, 0.65
                    bottom, height = 0.1, 0.65
                    bottom_h = left + width + 0.02
                    rect_scatter: Optional[Tuple[float, float, float, float]] = (
                        left,
                        bottom,
                        width,
                        height,
                    )
                    rect_histx = (left, bottom_h, width, 0.2)
                    # plt.set_facecolor('white')
                    nullfmt = NullFormatter()

                    axScatter = plt.axes(rect_scatter)
                    axScatter.set_facecolor("white")
                    p_array = np.array(p_values)
                    log2fc_array = np.array(log2fc_values)
                    log2fc_percentile = np.percentile(log2fc_array, 95)
                    ooFive = 0.05
                    ooOne = 0.01

                    # distribution
                    axHistx = plt.axes(rect_histx)
                    axHistx.set_facecolor("white")
                    axHistx.xaxis.set_major_formatter(nullfmt)
                    axHistx.yaxis.set_major_formatter(nullfmt)
                    binwidth = 0.05
                    xymax = np.max(
                        [np.max(np.fabs(log2fc_array)), np.max(np.fabs(p_values))]
                    )
                    lim = (int(xymax / binwidth) + 1) * binwidth
                    bins = np.arange(-lim, lim + binwidth, binwidth)
                    axHistx.hist(
                        log2fc_array,
                        bins=bins,
                        histtype="stepfilled",
                        color="grey",
                        align="mid",
                    )
                    # plot h-lines
                    # ax.axhline(y=ooFive, linewidth=2, color='orange', linestyle="--")
                    axScatter.axhline(
                        y=ooFive, linewidth=2, color="orange", linestyle="--"
                    )
                    ooFive_artist = Line2D(
                        (0, 1), (0, 0), color="orange", linestyle="--"
                    )
                    ooFive_label = f"p-value = {ooFive}"
                    ooOne_label = f"p-value = {ooOne}"
                    # ax.axhline(y=ooOne, linewidth=2, color='red', linestyle="--")
                    axScatter.axhline(y=ooOne, linewidth=2, color="red", linestyle="--")
                    ooOne_artist = Line2D((0, 1), (0, 0), color="red", linestyle="--")
                    # bonferroni
                    # ooFive_corrected = 0.05 / pair_data_count
                    # ooOne_corrected = 0.01 / pair_data_count
                    # ax.axhline(y=ooFive_corrected, linewidth=2, color='grey', linestyle="--")
                    # ooFive_corrected_artist = Line2D((0, 1), (0, 0), color='grey', linestyle='--')
                    # ax.axhline(y=ooOne_corrected, linewidth=2, color='black', linestyle="--")
                    # ooOne_corrected_artist = Line2D((0, 1), (0, 0), color='black', linestyle='--')

                    # plot v-lines
                    # ax.axvline(x=1.0, linewidth=2, color='purple', linestyle="--")
                    # ax.axvline(x=log2fc_percentile, linewidth=2, color='pink', linestyle="--")
                    axScatter.axvline(
                        x=1.0, linewidth=2, color="purple", linestyle="--"
                    )
                    axScatter.axvline(
                        x=log2fc_percentile, linewidth=2, color="blue", linestyle="--"
                    )
                    v1_label = "|log2FC| = 1"
                    v1_artist = Line2D((0, 1), (0, 0), color="purple", linestyle="--")
                    nine_five_percentile_label = "|log2FC-95%%ile| = %s" % (
                        "{0:.2f}".format(log2fc_percentile)
                    )
                    nine_five_percentile_artist = Line2D(
                        (0, 1), (0, 0), color="blue", linestyle="--"
                    )
                    # ax.axvline(x=-1.0, linewidth=2, color='purple', linestyle="--")
                    # ax.axvline(x=-log2fc_percentile, linewidth=2, color='pink', linestyle="--")
                    axScatter.axvline(
                        x=-1.0, linewidth=2, color="purple", linestyle="--"
                    )
                    axScatter.axvline(
                        x=-log2fc_percentile, linewidth=2, color="blue", linestyle="--"
                    )
                    # plot dots
                    # ax.scatter(log2fc_array, p_array, alpha=0.8, edgecolors='none', s=25, c='grey')
                    axScatter.scatter(
                        log2fc_array,
                        p_array,
                        alpha=0.8,
                        edgecolors="none",
                        s=25,
                        c="grey",
                    )
                    # Create legend from custom artist/label lists
                    # legend = ax.legend([ooFive_artist, ooOne_artist, ooFive_corrected_artist, ooOne_corrected_artist],
                    #          [ooFive, ooOne, "%s (0.05 corrected)" % '%.2E' % Decimal(ooFive_corrected), "%s (0.01 corrected)" % '%.2E' % Decimal(ooOne_corrected)],
                    #          fontsize=self.inputData.fontsize, frameon=True)
                    # legend = ax.legend([ooFive_artist, ooOne_artist, v1_artist, nine_five_percentile_artist],
                    legend = axScatter.legend(
                        [
                            ooFive_artist,
                            ooOne_artist,
                            v1_artist,
                            nine_five_percentile_artist,
                        ],
                        [
                            ooFive_label,
                            ooOne_label,
                            v1_label,
                            nine_five_percentile_label,
                        ],
                        fontsize=self.inputData.fontsize,
                        frameon=True,
                    )
                    legend.get_frame().set_facecolor("white")
                    if abs(np.min(log2fc_array)) < abs(np.max(log2fc_array)):
                        x_min = 0.0 - abs(np.max(log2fc_array))
                        x_max = 0.0 + abs(np.max(log2fc_array))
                    else:
                        x_min = 0.0 - abs(np.min(log2fc_array))
                        x_max = 0.0 + abs(np.min(log2fc_array))
                    # ax.set_xlim(x_min - 1, x_max + 1)
                    # ax.grid(True, linewidth=1, which="major", color="lightgrey")
                    # ax.grid(True, linewidth=0.5, which="minor", color="lightgrey")
                    # ax.set_ylim(np.min(p_array) * 0.1, 1.1)
                    # ax.set_xlabel("log2(mean(%s)/mean(%s))" % (x_label, y_label), fontsize=self.inputData.fontsize)
                    # ax.set_ylabel("p-value", fontsize=self.inputData.fontsize)
                    axScatter.set_xlim(x_min - 1, x_max + 1)
                    axScatter.grid(True, linewidth=1, which="major", color="lightgrey")
                    axScatter.grid(
                        True, linewidth=0.5, which="minor", color="lightgrey"
                    )
                    axScatter.set_ylim(1.1, np.min(p_array) * 0.1)
                    axScatter.set_xlabel(
                        f"log2(mean({x_label})/mean({y_label}))",
                        fontsize=self.inputData.fontsize,
                    )
                    axScatter.set_ylabel("p-value", fontsize=self.inputData.fontsize)
                    # ax.set_yscale('log')
                    axScatter.set_yscale("log")
                    axHistx.set_xlim(axScatter.get_xlim())
                    logger.info(f"[STATUS] - Plotting {pairwise_representation_test_f}")
                    # plt.gca().invert_yaxis()
                    plt.savefig(
                        pairwise_representation_test_f,
                        format=self.inputData.plot_format,
                    )
                    plt.close()

    def write_cluster_metrics(self) -> None:
        """
        Generate cluster metrics for the given data

        Parameters:
        - None

        Returns:
        - None: Metrics are saved as files.
        """
        cafe_f = os.path.join(self.dirs["main"], "cluster_counts_by_taxon.txt")
        cafe_output = []
        cafe_output.append(self.get_header_line("cafe", "TAXON"))

        cluster_metrics_domains_f = os.path.join(
            self.dirs["main"], "cluster_metrics_domains.txt"
        )
        cluster_metrics_domains_output = []
        cluster_metrics_domains_output.append(
            self.get_header_line("cluster_metrics_domains", "TAXON")
        )

        cluster_metrics_domains_detailed_output_by_domain_source = {}
        cluster_metrics_domains_detailed_f_by_domain_source = {}
        for domain_source in self.clusterCollection.domain_sources:
            cluster_metrics_domains_detailed_output_by_domain_source[domain_source] = []
            cluster_metrics_domains_detailed_output_by_domain_source[
                domain_source
            ].append(self.get_header_line("cluster_metrics_domains_detailed", "TAXON"))
            cluster_metrics_domains_detailed_f_by_domain_source[domain_source] = (
                os.path.join(
                    self.dirs["main"], f"cluster_domain_annotation.{domain_source}.txt"
                )
            )

        for attribute in self.aloCollection.attributes:

            attribute_metrics_f = os.path.join(
                self.dirs[attribute], f"{attribute}.attribute_metrics.txt"
            )
            attribute_metrics_output = []
            attribute_metrics_output.append(
                self.get_header_line("attribute_metrics", attribute)
            )

            pairwise_representation_test_f = os.path.join(
                self.dirs[attribute], f"{attribute}.pairwise_representation_test.txt"
            )
            pairwise_representation_test_output = []
            pairwise_representation_test_output.append(
                self.get_header_line("pairwise_representation_test", attribute)
            )

            pairwise_representation_test_by_pair_by_attribute = {}

            ###########################
            # cluster_metrics
            ###########################

            cluster_metrics_f = os.path.join(
                self.dirs[attribute], f"{attribute}.cluster_summary.txt"
            )
            cluster_metrics_output = []
            cluster_metrics_output.append(
                self.get_header_line("cluster_metrics", attribute)
            )

            levels = sorted(
                [x for x in self.aloCollection.ALO_by_level_by_attribute[attribute]]
            )
            levels_seen = set()

            for level in levels:
                ALO = self.aloCollection.ALO_by_level_by_attribute[attribute][level]
                if ALO:
                    ###########################
                    # attribute_metrics
                    ###########################

                    attribute_metrics_output.append(self.get_attribute_metrics(ALO))

                ###########################
                # cluster_metrics_ALO : setup
                ###########################

                cluster_metrics_ALO_f = os.path.join(
                    self.dirs[attribute], f"{attribute}.{level}.cluster_metrics.txt"
                )
                cluster_metrics_ALO_output = []
                cluster_metrics_ALO_output.append(
                    self.get_header_line("cluster_metrics_ALO", attribute)
                )

                background_representation_test_by_pair_by_attribute = {}

                ###########################
                # cluster_1to1s
                ###########################

                cluster_1to1_ALO_f = os.path.join(
                    self.dirs[attribute], f"{attribute}.{level}.cluster_1to1s.txt"
                )
                cluster_1to1_ALO_output = []
                cluster_1to1_ALO_output.append(
                    self.get_header_line("cluster_1to1s_ALO", attribute)
                )
                if not attribute == "TAXON":
                    if ALO:
                        for (
                            cluster_type
                        ) in ALO.clusters_by_cluster_cardinality_by_cluster_type:
                            for (
                                cluster_cardinality
                            ) in ALO.clusters_by_cluster_cardinality_by_cluster_type[
                                cluster_type
                            ]:
                                for (
                                    cluster_id
                                ) in ALO.clusters_by_cluster_cardinality_by_cluster_type[
                                    cluster_type
                                ][
                                    cluster_cardinality
                                ]:
                                    cluster_1to1_ALO_line = []
                                    cluster_1to1_ALO_line.append(cluster_id)
                                    cluster_1to1_ALO_line.append(cluster_type)
                                    cluster_1to1_ALO_line.append(cluster_cardinality)
                                    cluster_1to1_ALO_line.append(
                                        self.clusterCollection.cluster_list_by_cluster_id[
                                            cluster_id
                                        ].proteome_count
                                    )
                                    cluster_1to1_ALO_line.append(
                                        "{0:.2f}".format(
                                            len(
                                                [
                                                    protein_count
                                                    for proteome_id, protein_count in list(
                                                        self.clusterCollection.cluster_list_by_cluster_id[
                                                            cluster_id
                                                        ].protein_count_by_proteome_id.items()
                                                    )
                                                    if protein_count
                                                    == self.inputData.fuzzy_count
                                                ]
                                            )
                                            / self.clusterCollection.cluster_list_by_cluster_id[
                                                cluster_id
                                            ].proteome_count
                                        )
                                    )
                                    cluster_1to1_ALO_output.append(
                                        "\t".join(
                                            [
                                                str(field)
                                                for field in cluster_1to1_ALO_line
                                            ]
                                        )
                                    )

                for cluster in self.clusterCollection.cluster_list:

                    ###########################
                    # cluster_metrics (only done once for each attribute)
                    ###########################

                    if not levels_seen:
                        cluster_metrics_line = []
                        cluster_metrics_line.append(cluster.cluster_id)
                        cluster_metrics_line.append(cluster.protein_count)
                        cluster_metrics_line.append(cluster.protein_median)
                        cluster_metrics_line.append(cluster.proteome_count)
                        cluster_metrics_line.append(attribute)
                        cluster_metrics_line.append(
                            cluster.cluster_type_by_attribute[attribute]
                        )
                        if (
                            self.clusterCollection.fastas_parsed
                            and cluster.protein_length_stats
                        ):
                            cluster_metrics_line.append(
                                cluster.protein_length_stats["mean"]
                            )
                            cluster_metrics_line.append(
                                cluster.protein_length_stats["sd"]
                            )
                        else:
                            cluster_metrics_line.append("N/A")
                            cluster_metrics_line.append("N/A")
                        for _level in levels:
                            cluster_metrics_line.append(
                                sum(
                                    cluster.protein_counts_of_proteomes_by_level_by_attribute[
                                        attribute
                                    ][
                                        _level
                                    ]
                                )
                            )
                        if not attribute == "TAXON":
                            for _level in levels:
                                cluster_metrics_line.append(
                                    median(
                                        cluster.protein_counts_of_proteomes_by_level_by_attribute[
                                            attribute
                                        ][
                                            _level
                                        ]
                                    )
                                )
                            for _level in levels:
                                cluster_metrics_line.append(
                                    "{0:.2f}".format(
                                        cluster.proteome_coverage_by_level_by_attribute[
                                            attribute
                                        ][_level]
                                    )
                                )
                        cluster_metrics_output.append(
                            "\t".join([str(field) for field in cluster_metrics_line])
                        )

                    ###########################
                    # cafe (only done for attribute "TAXON")
                    ###########################

                    if not levels_seen and attribute == "TAXON":
                        cafe_line = []
                        # cafe_line.append("None")
                        cafe_line.append(str(cluster.cluster_id))
                        for _level in levels:
                            cafe_line.append(
                                sum(
                                    cluster.protein_counts_of_proteomes_by_level_by_attribute[
                                        attribute
                                    ][
                                        _level
                                    ]
                                )
                            )
                        cafe_output.append(
                            "\t".join([str(field) for field in cafe_line])
                        )

                    ###########################
                    # cluster_metrics_domains (only done for attribute "TAXON")
                    # - now different:
                    # - has line for each domain_id for each domain_source
                    ###########################

                    if not levels_seen and attribute == "TAXON":
                        if self.clusterCollection.functional_annotation_parsed:
                            # cluster_metrics_domain_line
                            cluster_metrics_domains_line = []
                            cluster_metrics_domains_line.append(cluster.cluster_id)
                            cluster_metrics_domains_line.append(cluster.protein_count)
                            cluster_metrics_domains_line.append(cluster.proteome_count)
                            if (
                                self.clusterCollection.fastas_parsed
                                and cluster.protein_length_stats
                            ):
                                cluster_metrics_domains_line.append(
                                    cluster.protein_length_stats["mean"]
                                )
                                cluster_metrics_domains_line.append(
                                    cluster.protein_length_stats["sd"]
                                )
                            else:
                                cluster_metrics_domains_line.append("N/A")
                                cluster_metrics_domains_line.append("N/A")
                            if "SignalP_EUK" in self.clusterCollection.domain_sources:
                                cluster_metrics_domains_line.append(
                                    "{0:.2f}".format(cluster.secreted_cluster_coverage)
                                )
                            else:
                                cluster_metrics_domains_line.append("N/A")
                            for domain_source in self.clusterCollection.domain_sources:
                                # cluster_metrics_domains
                                if (
                                    domain_source
                                    in cluster.domain_counter_by_domain_source
                                ):
                                    sorted_counts = sorted(
                                        [
                                            f"{domain_id}:{count}"
                                            for domain_id, count in cluster.domain_counter_by_domain_source[
                                                domain_source
                                            ].most_common()
                                        ],
                                        key=lambda x: (
                                            x.split(":")[-1],
                                            x.split(":")[-2],
                                        ),
                                    )
                                    sorted_counts_str = ";".join(sorted_counts)
                                    cluster_metrics_domains_line.append(
                                        sorted_counts_str
                                    )
                                    cluster_metrics_domains_line.append(
                                        "{0:.3f}".format(
                                            cluster.domain_entropy_by_domain_source[
                                                domain_source
                                            ]
                                        )
                                    )
                                else:
                                    cluster_metrics_domains_line.append("N/A")
                                    cluster_metrics_domains_line.append("N/A")
                            cluster_metrics_domains_output.append(
                                "\t".join(
                                    [
                                        str(field)
                                        for field in cluster_metrics_domains_line
                                    ]
                                )
                            )
                            for (
                                domain_source
                            ) in cluster.domain_counter_by_domain_source:
                                for (
                                    domain_id,
                                    count,
                                ) in cluster.domain_counter_by_domain_source[
                                    domain_source
                                ].most_common():
                                    cluster_metrics_domains_detailed_output_line = []
                                    cluster_metrics_domains_detailed_output_line.append(
                                        cluster.cluster_id
                                    )
                                    cluster_metrics_domains_detailed_output_line.append(
                                        domain_source
                                    )
                                    cluster_metrics_domains_detailed_output_line.append(
                                        domain_id
                                    )
                                    if domain_source == "SignalP_EUK":
                                        cluster_metrics_domains_detailed_output_line.append(
                                            domain_id
                                        )
                                    else:
                                        if (
                                            domain_source
                                            in self.proteinCollection.domain_desc_by_id_by_source
                                        ):
                                            cluster_metrics_domains_detailed_output_line.append(
                                                self.proteinCollection.domain_desc_by_id_by_source[
                                                    domain_source
                                                ].get(
                                                    domain_id, "N/A"
                                                )
                                            )
                                        else:
                                            cluster_metrics_domains_detailed_output_line.append(
                                                "N/A"
                                            )
                                    cluster_metrics_domains_detailed_output_line.append(
                                        cluster.protein_count
                                    )
                                    protein_with_domain_count_by_proteome_id = {}
                                    proteome_count_with_domain = 0
                                    protein_without_domain_count_by_proteome_id = {}
                                    for proteome_id, protein_ids in list(
                                        cluster.protein_ids_by_proteome_id.items()
                                    ):
                                        proteome_seen = False
                                        for protein_id in protein_ids:
                                            if (
                                                domain_source
                                                in self.proteinCollection.proteins_by_protein_id[
                                                    protein_id
                                                ].domain_counter_by_domain_source
                                                and domain_id
                                                in self.proteinCollection.proteins_by_protein_id[
                                                    protein_id
                                                ].domain_counter_by_domain_source[
                                                    domain_source
                                                ]
                                            ):
                                                protein_with_domain_count_by_proteome_id[
                                                    proteome_id
                                                ] = (
                                                    protein_with_domain_count_by_proteome_id.get(
                                                        proteome_id, 0
                                                    )
                                                    + 1
                                                )
                                                if not proteome_seen:
                                                    proteome_count_with_domain += 1
                                                    proteome_seen = True
                                            else:
                                                protein_without_domain_count_by_proteome_id[
                                                    proteome_id
                                                ] = (
                                                    protein_without_domain_count_by_proteome_id.get(
                                                        proteome_id, 0
                                                    )
                                                    + 1
                                                )
                                    proteomes_with_domain_count_string = ",".join(
                                        sorted(
                                            [
                                                f"{proteome_id}:{count}/{len(cluster.protein_ids_by_proteome_id[proteome_id])}"
                                                for proteome_id, count in list(
                                                    protein_with_domain_count_by_proteome_id.items()
                                                )
                                            ]
                                        )
                                    )
                                    proteomes_without_domain_count_string = ",".join(
                                        sorted(
                                            [
                                                f"{proteome_id}:{count}/{len(cluster.protein_ids_by_proteome_id[proteome_id])}"
                                                for proteome_id, count in list(
                                                    protein_without_domain_count_by_proteome_id.items()
                                                )
                                            ]
                                        )
                                    )
                                    cluster_metrics_domains_detailed_output_line.append(
                                        sum(
                                            protein_with_domain_count_by_proteome_id.values()
                                        )
                                    )
                                    cluster_metrics_domains_detailed_output_line.append(
                                        "{0:.3f}".format(
                                            proteome_count_with_domain
                                            / cluster.proteome_count
                                        )
                                    )
                                    if proteomes_with_domain_count_string:
                                        cluster_metrics_domains_detailed_output_line.append(
                                            proteomes_with_domain_count_string
                                        )
                                    else:
                                        cluster_metrics_domains_detailed_output_line.append(
                                            "N/A"
                                        )
                                    if proteomes_without_domain_count_string:
                                        cluster_metrics_domains_detailed_output_line.append(
                                            proteomes_without_domain_count_string
                                        )
                                    else:
                                        cluster_metrics_domains_detailed_output_line.append(
                                            "N/A"
                                        )
                                    cluster_metrics_domains_detailed_output_by_domain_source[
                                        domain_source
                                    ].append(
                                        "\t".join(
                                            [
                                                str(field)
                                                for field in cluster_metrics_domains_detailed_output_line
                                            ]
                                        )
                                    )

                    ###########################
                    # cluster_metrics_ALO : populate
                    ###########################

                    cluster_metrics_ALO_line = []
                    cluster_metrics_ALO_line.append(cluster.cluster_id)
                    if ALO:
                        cluster_metrics_ALO_line.append(
                            ALO.cluster_status_by_cluster_id[cluster.cluster_id]
                        )
                        cluster_metrics_ALO_line.append(
                            ALO.cluster_type_by_cluster_id[cluster.cluster_id]
                        )
                    cluster_metrics_ALO_line.append(cluster.protein_count)
                    cluster_metrics_ALO_line.append(cluster.proteome_count)
                    cluster_metrics_ALO_line.append(
                        sum(
                            cluster.protein_counts_of_proteomes_by_level_by_attribute[
                                attribute
                            ][level]
                        )
                    )
                    if (
                        ALO
                        and ALO.cluster_mean_ALO_count_by_cluster_id[cluster.cluster_id]
                    ):
                        cluster_metrics_ALO_line.append(
                            ALO.cluster_mean_ALO_count_by_cluster_id[cluster.cluster_id]
                        )
                    else:
                        cluster_metrics_ALO_line.append("N/A")
                    if (
                        ALO
                        and ALO.cluster_mean_non_ALO_count_by_cluster_id[
                            cluster.cluster_id
                        ]
                    ):
                        cluster_metrics_ALO_line.append(
                            ALO.cluster_mean_non_ALO_count_by_cluster_id[
                                cluster.cluster_id
                            ]
                        )
                    else:
                        cluster_metrics_ALO_line.append("N/A")
                    if (
                        ALO
                        and ALO.cluster_type_by_cluster_id[cluster.cluster_id]
                        == "shared"
                    ):
                        if ALO.cluster_mwu_log2_mean_by_cluster_id[cluster.cluster_id]:
                            background_pair = (level, "background")
                            if (
                                attribute
                                not in background_representation_test_by_pair_by_attribute
                            ):
                                background_representation_test_by_pair_by_attribute[
                                    attribute
                                ] = {}
                            if (
                                background_pair
                                not in background_representation_test_by_pair_by_attribute[
                                    attribute
                                ]
                            ):
                                background_representation_test_by_pair_by_attribute[
                                    attribute
                                ][background_pair] = []
                            background_representation_test = []
                            background_representation_test.append(cluster.cluster_id)
                            background_representation_test.append(level)
                            background_representation_test.append("background")
                            background_representation_test.append(
                                ALO.cluster_mean_ALO_count_by_cluster_id[
                                    cluster.cluster_id
                                ]
                            )
                            background_representation_test.append(
                                ALO.cluster_mean_non_ALO_count_by_cluster_id[
                                    cluster.cluster_id
                                ]
                            )
                            background_representation_test.append(
                                ALO.cluster_mwu_log2_mean_by_cluster_id[
                                    cluster.cluster_id
                                ]
                            )
                            background_representation_test.append(
                                ALO.cluster_mwu_pvalue_by_cluster_id[cluster.cluster_id]
                            )
                            background_representation_test_by_pair_by_attribute[
                                attribute
                            ][background_pair].append(background_representation_test)

                            if (
                                ALO.cluster_mwu_log2_mean_by_cluster_id[
                                    cluster.cluster_id
                                ]
                                > 0
                            ):
                                cluster_metrics_ALO_line.append("enriched")
                            elif (
                                ALO.cluster_mwu_log2_mean_by_cluster_id[
                                    cluster.cluster_id
                                ]
                                < 0
                            ):
                                cluster_metrics_ALO_line.append("depleted")
                            else:
                                cluster_metrics_ALO_line.append("equal")
                            cluster_metrics_ALO_line.append(
                                ALO.cluster_mwu_log2_mean_by_cluster_id[
                                    cluster.cluster_id
                                ]
                            )
                            cluster_metrics_ALO_line.append(
                                ALO.cluster_mwu_pvalue_by_cluster_id[cluster.cluster_id]
                            )
                        else:
                            cluster_metrics_ALO_line.append("N/A")
                            cluster_metrics_ALO_line.append("N/A")
                            cluster_metrics_ALO_line.append("N/A")
                    else:
                        cluster_metrics_ALO_line.append("N/A")
                        cluster_metrics_ALO_line.append("N/A")
                        cluster_metrics_ALO_line.append("N/A")
                    cluster_metrics_ALO_line.append(
                        "{0:.2f}".format(
                            cluster.proteome_coverage_by_level_by_attribute[attribute][
                                level
                            ]
                        )
                    )
                    ALO_proteomes_present = cluster.proteome_ids.intersection(
                        ALO.proteomes if ALO else set("")
                    )
                    non_ALO_proteomes_present = cluster.proteome_ids.difference(
                        ALO.proteomes if ALO else set("")
                    )
                    cluster_metrics_ALO_line.append(len(ALO_proteomes_present))
                    cluster_metrics_ALO_line.append(len(non_ALO_proteomes_present))
                    if ALO_proteomes_present:
                        cluster_metrics_ALO_line.append(
                            ",".join(sorted(list(ALO_proteomes_present)))
                        )
                    else:
                        cluster_metrics_ALO_line.append("N/A")
                    if non_ALO_proteomes_present:
                        cluster_metrics_ALO_line.append(
                            ",".join(sorted(list(non_ALO_proteomes_present)))
                        )
                    else:
                        cluster_metrics_ALO_line.append("N/A")
                    # if cluster.go_terms:
                    #    cluster_metrics_ALO_line.append(";".join(sorted(list(cluster.go_terms))))
                    # else:
                    #    cluster_metrics_ALO_line.append("N/A")
                    # for domain_source in self.clusterCollection.domain_sources:
                    #    if domain_source in cluster.domain_counter_by_domain_source:
                    #        cluster_metrics_ALO_line.append(";".join(["%s:%s" % (domain, count) for domain, count in cluster.domain_counter_by_domain_source[domain_source].most_common()]))
                    #    else:
                    #        cluster_metrics_ALO_line.append("N/A")
                    cluster_metrics_ALO_output.append(
                        "\t".join([str(field) for field in cluster_metrics_ALO_line])
                    )

                    if (
                        len(levels) > 1
                        and len(ALO_proteomes_present) >= self.inputData.min_proteomes
                    ):
                        for result in self.pairwise_representation_test(
                            cluster, attribute, level, levels_seen, levels
                        ):
                            # [cluster.cluster_id, level, other_level, mean_level, mean_other_level, log2fc_mean, pvalue]
                            if (
                                attribute
                                not in pairwise_representation_test_by_pair_by_attribute
                            ):
                                pairwise_representation_test_by_pair_by_attribute[
                                    attribute
                                ] = {}
                            pair = (result[1], result[2])
                            if (
                                pair
                                not in pairwise_representation_test_by_pair_by_attribute[
                                    attribute
                                ]
                            ):
                                pairwise_representation_test_by_pair_by_attribute[
                                    attribute
                                ][pair] = []
                            pairwise_representation_test_by_pair_by_attribute[
                                attribute
                            ][pair].append(result)

                            pairwise_representation_test_line = []
                            pairwise_representation_test_line.append(result[0])
                            pairwise_representation_test_line.append(result[1])
                            pairwise_representation_test_line.append(result[3])
                            pairwise_representation_test_line.append(result[2])
                            pairwise_representation_test_line.append(result[4])
                            pairwise_representation_test_line.append(result[5])
                            pairwise_representation_test_line.append(result[6])
                            # if cluster.go_terms:
                            #    pairwise_representation_test_line.append(";".join(sorted(list(cluster.go_terms))))
                            # else:
                            #    pairwise_representation_test_line.append("N/A")
                            # for domain_source in clusterCollection.domain_sources:
                            #    if domain_source in cluster.domain_counter_by_domain_source:
                            #        pairwise_representation_test_line.append(";".join(["%s:%s" % (domain, count) for domain, count in cluster.domain_counter_by_domain_source[domain_source].most_common()]))
                            #    else:
                            #        pairwise_representation_test_line.append("N/A")
                            pairwise_representation_test_output.append(
                                "\t".join(
                                    [
                                        str(field)
                                        for field in pairwise_representation_test_line
                                    ]
                                )
                            )

                levels_seen.add(level)
                # END of cluster loop

                if len(cafe_output) > 1:
                    with open(cafe_f, "w") as cafe_fh:
                        logger.info(f"[STATUS] - Writing {cafe_f}")
                        cafe_fh.write("\n".join(cafe_output) + "\n")
                    cafe_output = []
                if len(cluster_metrics_output) > 1:
                    with open(cluster_metrics_f, "w") as cluster_metrics_fh:
                        logger.info(f"[STATUS] - Writing {cluster_metrics_f}")
                        cluster_metrics_fh.write(
                            "\n".join(cluster_metrics_output) + "\n"
                        )
                    cluster_metrics_output = []
                if len(cluster_metrics_domains_output) > 1:
                    with open(
                        cluster_metrics_domains_f, "w"
                    ) as cluster_metrics_domains_fh:
                        logger.info(f"[STATUS] - Writing {cluster_metrics_domains_f}")
                        cluster_metrics_domains_fh.write(
                            "\n".join(cluster_metrics_domains_output) + "\n"
                        )
                    cluster_metrics_domains_output = []
                for (
                    domain_source
                ) in cluster_metrics_domains_detailed_output_by_domain_source:
                    if (
                        len(
                            cluster_metrics_domains_detailed_output_by_domain_source[
                                domain_source
                            ]
                        )
                        > 1
                    ):
                        cluster_metrics_domains_detailed_f = (
                            cluster_metrics_domains_detailed_f_by_domain_source[
                                domain_source
                            ]
                        )
                        with open(
                            cluster_metrics_domains_detailed_f, "w"
                        ) as cluster_metrics_domains_detailed_fh:
                            logger.info(
                                f"[STATUS] - Writing {cluster_metrics_domains_detailed_f}"
                            )
                            cluster_metrics_domains_detailed_fh.write(
                                "\n".join(
                                    cluster_metrics_domains_detailed_output_by_domain_source[
                                        domain_source
                                    ]
                                )
                                + "\n"
                            )
                        cluster_metrics_domains_detailed_output_by_domain_source[
                            domain_source
                        ] = []
                if len(cluster_metrics_ALO_output) > 1:
                    with open(cluster_metrics_ALO_f, "w") as cluster_metrics_ALO_fh:
                        logger.info(f"[STATUS] - Writing {cluster_metrics_ALO_f}")
                        cluster_metrics_ALO_fh.write(
                            "\n".join(cluster_metrics_ALO_output) + "\n"
                        )
                    cluster_metrics_ALO_output = []
                if len(cluster_1to1_ALO_output) > 1:
                    with open(cluster_1to1_ALO_f, "w") as cluster_1to1_ALO_fh:
                        logger.info(f"[STATUS] - Writing {cluster_1to1_ALO_f}")
                        cluster_1to1_ALO_fh.write(
                            "\n".join(cluster_1to1_ALO_output) + "\n"
                        )
                    cluster_1to1_ALO_output = []
                if background_representation_test_by_pair_by_attribute:
                    self.plot_count_comparisons_volcano(
                        background_representation_test_by_pair_by_attribute
                    )

            if len(attribute_metrics_output) > 1:
                with open(attribute_metrics_f, "w") as attribute_metrics_fh:
                    logger.info(f"[STATUS] - Writing {attribute_metrics_f}")
                    attribute_metrics_fh.write(
                        "\n".join(attribute_metrics_output) + "\n"
                    )
            if len(pairwise_representation_test_output) > 1:
                with open(
                    pairwise_representation_test_f, "w"
                ) as pairwise_representation_test_fh:
                    logger.info(f"[STATUS] - Writing {pairwise_representation_test_f}")
                    pairwise_representation_test_fh.write(
                        "\n".join(pairwise_representation_test_output) + "\n"
                    )
            if pairwise_representation_test_by_pair_by_attribute:
                self.plot_count_comparisons_volcano(
                    pairwise_representation_test_by_pair_by_attribute
                )

    def write_output(self) -> None:
        """
        generates the results and stores them into files
        """
        self.plot_cluster_sizes()
        self.write_cluster_metrics()
