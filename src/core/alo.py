from typing import Dict, List, Literal, Optional, Set


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

        self.rarefaction_data = {}  # repetition : number of clusters
