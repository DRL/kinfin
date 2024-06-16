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

    def add_annotation_to_protein(
        self,
        domain_protein_id: str,
        domain_counter_by_domain_source: Dict[str, Counter],
        go_terms: List[str],
    ):
        """
        Updates a protein object with domain counters and GO terms.

        Args:
        - domain_protein_id (str): Identifier of the protein to annotate.
        - domain_counter_by_domain_source (Dict[str, Counter]): Domain sources mapped to counters of domains.
        - go_terms (List[str]): Gene Ontology (GO) terms associated with the protein.

        This method sets domain counters, assigns GO terms, and checks if the protein is secreted
        based on domain information ('SignalP_EUK' source).

        Note: If 'SignalP_EUK' indicates 'SignalP-noTM', sets protein.secreted = True.
        """
        protein: Optional[Protein] = self.proteins_by_protein_id.get(
            domain_protein_id, None
        )
        if protein is not None:
            protein.domain_counter_by_domain_source = domain_counter_by_domain_source
            signalp_notm = protein.domain_counter_by_domain_source.get(
                "SignalP_EUK", None
            )
            if signalp_notm and "SignalP-noTM" in signalp_notm:
                protein.secreted = True
            protein.go_terms = go_terms
