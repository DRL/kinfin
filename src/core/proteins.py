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
