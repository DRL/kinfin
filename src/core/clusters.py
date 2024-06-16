from collections import Counter
from typing import Dict, FrozenSet, List

from core.proteins import ProteinCollection


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
