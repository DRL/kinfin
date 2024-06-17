from collections import Counter
from math import log
from typing import DefaultDict, Dict, FrozenSet, List, Literal, Optional, Set

from core.logic import compute_protein_ids_by_proteome
from core.proteins import ProteinCollection
from core.utils import mean, median, sd


class Cluster:
    def __init__(
        self,
        cluster_id: str,
        protein_ids: List[str],
        proteinCollection: ProteinCollection,
    ) -> None:
        self.cluster_id: str = cluster_id
        self.protein_ids = set(protein_ids)
        self.protein_count: int = len(protein_ids)
        try:

            self.proteomes_by_protein_id: Dict[str, str] = {
                id: proteinCollection.proteins_by_protein_id[id].proteome_id
                for id in protein_ids
            }
        except KeyError as e:
            error_msg = f"[ERROR] - Protein {e.args[0]} in clustering belongs to proteomes that are not present in the config-file. Please add those proteomes or recluster by omitting these proteomes."
            raise KeyError(error_msg)

        self.proteome_ids_list: List[str] = list(self.proteomes_by_protein_id.values())
        self.protein_count_by_proteome_id: Counter[str] = Counter(self.proteome_ids_list)  # fmt: skip
        self.proteome_ids: FrozenSet[str] = frozenset(self.proteome_ids_list)
        self.proteome_count: int = len(self.proteome_ids)
        self.singleton: bool = False if self.protein_count > 1 else True
        self.apomorphy: bool = False if self.proteome_count > 1 else True

        self.protein_ids_by_proteome_id: DefaultDict[str, Set[str]] = compute_protein_ids_by_proteome(self.proteomes_by_protein_id)  # fmt:skip
        self.protein_counts_of_proteomes_by_level_by_attribute:Dict[str, Dict[str, List[int]]] = {}  # fmt:skip
        self.proteome_coverage_by_level_by_attribute: Dict[str, Dict[str, float]] = {}
        self.implicit_protein_ids_by_proteome_id_by_level_by_attribute: Dict[str, Dict[str, Dict[str, List[str]]]] = {}  # fmt:skip
        self.cluster_type_by_attribute: Dict[
            str,
            Literal["singleton", "shared", "specific"],
        ] = {}
        self.protein_median: Optional[float] = None
        self.protein_length_stats: Optional[Dict[str,float]]= self.compute_protein_length_stats(proteinCollection, self.protein_ids)  # fmt:skip
        self.secreted_cluster_coverage: float = self.compute_secreted_cluster_coverage(proteinCollection, self.protein_ids, self.protein_count)  # fmt:skip
        self.domain_counter_by_domain_source: Dict[str, Counter[str]] = self.compute_domain_counter_by_domain_source(proteinCollection, self.protein_ids)  # fmt:skip
        self.domain_entropy_by_domain_source: Dict[str,float] = self.compute_domain_entropy_by_domain_source()  # fmt:skip

    def compute_protein_length_stats(
        self,
        proteinCollection: ProteinCollection,
        protein_ids: Set[str],
    ) -> Optional[Dict[str, float]]:
        """
        Computes statistics (mean, median, standard deviation) of protein lengths.

        Parameters:
        - proteinCollection: A ProteinCollection object containing protein data.
        - protein_ids: A set of protein IDs for which lengths are to be computed.

        Returns:
        - Optional[Dict[str, float]]: A dictionary containing 'mean', 'median', and 'sd'
          (standard deviation) of protein lengths, if all lengths are available and at least
          one protein ID is provided. Returns None if no valid protein lengths are found.
        """
        protein_lengths: List[int | None] = [
            proteinCollection.proteins_by_protein_id[protein_id].length
            for protein_id in protein_ids
        ]
        if all(protein_lengths):
            protein_length_stats: Dict[str, float] = {}
            protein_length_stats["mean"] = mean(protein_lengths)
            protein_length_stats["median"] = median(protein_lengths)
            protein_length_stats["sd"] = sd(protein_lengths)
            return protein_length_stats

    def compute_secreted_cluster_coverage(
        self,
        proteinCollection: ProteinCollection,
        protein_ids: Set[str],
        protein_count: int,
    ) -> float:
        """
        Computes the fraction of secreted proteins in a given set of protein IDs.

        Parameters:
        - proteinCollection: A ProteinCollection object containing protein data.
        - protein_ids: A set of protein IDs to compute secreted protein coverage.
        - protein_count: Total count of proteins in the cluster.

        Returns:
        - float: Fraction of secreted proteins in the provided set of protein IDs.
        """
        secreted = 0
        for protein_id in protein_ids:
            if proteinCollection.proteins_by_protein_id[protein_id].secreted:
                secreted += 1
        return secreted / protein_count

    def compute_domain_counter_by_domain_source(
        self,
        proteinCollection: ProteinCollection,
        protein_ids: Set[str],
    ) -> Dict[str, Counter[str]]:
        """
        Computes the aggregated domain counts by domain source for a set of protein IDs.

        Parameters:
        - proteinCollection: A ProteinCollection object containing protein data.
        - protein_ids: A set of protein IDs for which domain counts are computed.

        Returns:
        - Dict[str, Counter[str]]: A dictionary where keys are domain sources and values are
          Counters mapping domain IDs to their respective counts.
        """
        cluster_domain_counter_by_domain_source: Dict[str, Counter[str]] = {}
        for protein_id in protein_ids:
            protein_domain_counter_by_domain_source: Dict[str, Counter[str]] = (
                proteinCollection.proteins_by_protein_id[
                    protein_id
                ].domain_counter_by_domain_source
            )
            if protein_domain_counter_by_domain_source:
                for domain_source, protein_domain_counter in list(
                    protein_domain_counter_by_domain_source.items()
                ):
                    if domain_source not in cluster_domain_counter_by_domain_source:
                        cluster_domain_counter_by_domain_source[domain_source] = (
                            Counter()
                        )
                    cluster_domain_counter_by_domain_source[
                        domain_source
                    ] += protein_domain_counter
        return cluster_domain_counter_by_domain_source

    def compute_domain_entropy_by_domain_source(self) -> Dict[str, float]:
        """
        Computes entropy for domains grouped by different sources.

        Returns:
        - Dict[str, float]: Dictionary where keys are domain sources and values are computed entropy values.
        """
        self.domain_entropy_by_domain_source: Dict[str, float] = {}
        for domain_source, domain_counter in list(
            self.domain_counter_by_domain_source.items()
        ):
            total_count: int = len([domain for domain in domain_counter.elements()])
            domain_entropy: float = -sum(
                [
                    i / total_count * log(i / total_count, 2)
                    for i in list(domain_counter.values())
                ]
            )
            if str(domain_entropy) == "-0.0":
                self.domain_entropy_by_domain_source[domain_source] = 0.0
            else:
                self.domain_entropy_by_domain_source[domain_source] = domain_entropy
        return self.domain_entropy_by_domain_source


class ClusterCollection:
    def __init__(
        self,
        cluster_list: List[Cluster],
        inferred_singletons_count: int,
        functional_annotation_parsed: bool,
        fastas_parsed: bool,
        domain_sources: List[str],
    ):
        self.cluster_list: List[Cluster] = cluster_list
        self.cluster_list_by_cluster_id: Dict[str, Cluster] = {
            cluster.cluster_id: cluster for cluster in cluster_list
        }  # only for testing
        self.cluster_count: int = len(cluster_list)
        self.inferred_singletons_count: int = inferred_singletons_count
        self.functional_annotation_parsed: bool = functional_annotation_parsed
        self.fastas_parsed: bool = fastas_parsed
        # self.domain_sources = [domain_source for domain_source in domain_sources if not domain_source == "GO"]
        self.domain_sources: List[str] = domain_sources
