from collections import Counter
from typing import Dict, List, Optional


class Protein:
    def __init__(
        self,
        protein_id: str,
        proteome_id: str,
        species_id: str,
        sequence_id: str,
    ) -> None:

        self.protein_id = protein_id
        self.proteome_id = proteome_id
        self.species_id = species_id
        self.sequence_id = sequence_id
        self.length: Optional[int] = None
        self.clustered: bool = False
        self.secreted: bool = False
        self.domain_counter_by_domain_source: Dict[str, Counter[str]] = {}
        self.go_terms: List[str] = []

    def update_length(self, length: int) -> None:
        self.length = length


class ProteinCollection:
    def __init__(self, proteins_list: List[Protein]) -> None:
        self.proteins_list: List[Protein] = proteins_list
        self.proteins_by_protein_id: Dict[str, Protein] = {
            protein.protein_id: protein for protein in proteins_list
        }
        self.protein_count: int = len(proteins_list)
        self.domain_sources: List[str] = []
        self.fastas_parsed: bool = False
        self.functional_annotation_parsed: bool = False
        self.domain_desc_by_id_by_source: Dict[str, Dict[str,str]] = {}  # fmt: skip
