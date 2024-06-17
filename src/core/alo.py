from typing import Dict, List, Literal, Optional, Set

from core.clusters import Cluster


class AttributeLevel:
    """
    Definitions:
        'shared' : shared between one ALO and others
        'singleton' : cardinality of 1 ('specific', but separate)
        'specific' : only present within one ALO
    """

    def __init__(self, attribute: str, level: str, proteomes: Set[str]) -> None:
        self.attribute: str = attribute
        self.level: str = level
        self.proteomes: Set[str] = set(proteomes)
        self.proteomes_list: List[str] = list(proteomes)
        self.proteome_count: int = len(proteomes)

        self.cluster_ids_by_cluster_type_by_cluster_status: Dict[
            str, Dict[str, List[str]]
        ] = {
            # sums up to cluster_count
            "present": {"singleton": [], "specific": [], "shared": []},
            "absent": {"singleton": [], "specific": [], "shared": []},
        }

        self.protein_ids_by_cluster_type: Dict[str, List[str]] = {
            # list of lists
            "singleton": [],
            "specific": [],
            "shared": [],
        }

        self.protein_span_by_cluster_type: Dict[str, List[int | float]] = {
            "singleton": [],
            "specific": [],
            "shared": [],
        }

        self.clusters_by_cluster_cardinality_by_cluster_type: Dict[
            str, Dict[str, List[str]]
        ] = {
            "shared": {"true": [], "fuzzy": []},
            "specific": {"true": [], "fuzzy": []},
        }

        self.cluster_status_by_cluster_id: Dict[str, Literal["absent", "present"]] = {}
        self.cluster_type_by_cluster_id: Dict[str,Literal['singleton', 'shared', 'specific']] = {}  # fmt:skip

        self.cluster_mwu_pvalue_by_cluster_id = {}
        self.cluster_mwu_log2_mean_by_cluster_id = {}
        self.cluster_mean_ALO_count_by_cluster_id = {}
        self.cluster_mean_non_ALO_count_by_cluster_id = {}

        self.domain_counter_by_domain_source_by_cluster_type = None
        self.protein_with_domain_count_by_domain_source_by_cluster_type = None

        self.protein_length_stats_by_cluster_id: Dict[str, Dict[str, int | float]] = {}  # fmt:skip
        self.protein_count_by_cluster_id: Dict[str, int] = {}

    def add_cluster(
        self,
        cluster: Cluster,
        attribute_cluster_type: Literal["singleton", "shared", "specific"],
        ALO_cluster_status: Literal["absent", "present"],
        ALO_protein_length_stats: Dict[str, int | float],
        ALO_protein_ids_in_cluster: List[str],
        ALO_cluster_cardinality: Optional[str],
        mwu_pvalue: Optional[float],
        mwu_log2_mean: Optional[float],
        mean_ALO_count: Optional[float],
        mean_non_ALO_count: Optional[float],
    ) -> None:
        """
        Adds a cluster to various data structures maintained by the class.

        Args:
            cluster (Cluster): The cluster object to add.
            attribute_cluster_type (Literal["singleton", "shared", "specific"]):
                Type of the cluster as either 'singleton', 'shared', or 'specific'.
            ALO_cluster_status (Literal["absent", "present"]):
                Status of the cluster, either 'absent' or 'present'.
            ALO_protein_length_stats (Dict[str, int | float]):
                Length statistics of proteins in the cluster.
            ALO_protein_ids_in_cluster (List[str]):
                List of protein IDs present in the cluster.
            ALO_cluster_cardinality (Optional[str]):
                Cardinality of the cluster (if applicable).
            mwu_pvalue (Optional[float]):
                P-value from Mann-Whitney U test (if applicable).
            mwu_log2_mean (Optional[float]):
                Log2 transformed mean (if applicable).
            mean_ALO_count (Optional[float]):
                Mean count of ALO (if applicable).
            mean_non_ALO_count (Optional[float]):
                Mean count of non-ALO (if applicable).

        Returns:
            None
        """
        self.cluster_ids_by_cluster_type_by_cluster_status[ALO_cluster_status][
            attribute_cluster_type
        ].append(cluster.cluster_id)
        self.cluster_status_by_cluster_id[cluster.cluster_id] = ALO_cluster_status
        self.cluster_type_by_cluster_id[cluster.cluster_id] = attribute_cluster_type
        self.protein_length_stats_by_cluster_id[cluster.cluster_id] = ALO_protein_length_stats  # fmt:skip

        self.protein_count_by_cluster_id[cluster.cluster_id] = len(ALO_protein_ids_in_cluster)  # fmt:skip

        if ALO_cluster_status == "present":
            for ALO_protein_id in ALO_protein_ids_in_cluster:
                self.protein_ids_by_cluster_type[attribute_cluster_type].append(ALO_protein_id)  # fmt:skip
            self.protein_span_by_cluster_type[attribute_cluster_type].append(ALO_protein_length_stats["sum"])  # fmt:skip
            if not attribute_cluster_type == "singleton":
                if ALO_cluster_cardinality:
                    self.clusters_by_cluster_cardinality_by_cluster_type[
                        attribute_cluster_type
                    ][ALO_cluster_cardinality].append(cluster.cluster_id)

        self.cluster_mwu_pvalue_by_cluster_id[cluster.cluster_id] = mwu_pvalue
        self.cluster_mwu_log2_mean_by_cluster_id[cluster.cluster_id] = mwu_log2_mean
        self.cluster_mean_ALO_count_by_cluster_id[cluster.cluster_id] = mean_ALO_count
        self.cluster_mean_non_ALO_count_by_cluster_id[cluster.cluster_id] = (
            mean_non_ALO_count
        )
