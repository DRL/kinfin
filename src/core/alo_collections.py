from typing import Any, Dict, List, Optional, Set, Tuple

from ete3 import Tree

from core.config import ATTRIBUTE_RESERVED


class AloCollection:
    def __init__(
        self,
        proteomes: Set[str],
        attributes: List[str],
        proteome_id_by_species_id: Dict[str, str],
        level_by_attribute_by_proteome_id: Dict[str, Dict[str, str]],
        node_idx_by_proteome_ids: Optional[Dict[Any, Any]],
        tree_ete: Optional[Tree],
    ) -> None:
        self.proteomes = proteomes
        self.attributes_verbose = attributes
        self.attributes = [
            # list of attributes
            attribute
            for attribute in attributes
            if attribute not in ATTRIBUTE_RESERVED
        ]
        self.proteome_id_by_species_id = proteome_id_by_species_id
        self.level_by_attribute_by_proteome_id = level_by_attribute_by_proteome_id
        self.node_idx_by_proteome_ids = node_idx_by_proteome_ids
        self.tree_ete = tree_ete
        self.proteome_ids_by_level_by_attribute = self.compute_proteomes_by_level_by_attribute()  # fmt: skip

    def compute_proteomes_by_level_by_attribute(self) -> Dict[str, Dict[str, Set[str]]]:
        """
        Compute proteomes grouped by levels for each attribute.

        Args:
            attributes (List[str]): A list of strings representing attributes.
            level_by_attribute_by_proteome_id (Dict[str, Dict[str, str]]): A dictionary where keys
                are proteome IDs (strings), and values are dictionaries with keys representing
                attributes (strings) and values representing levels (strings).

        Returns:
            Dict[str, Dict[str, Set[str]]]: A dictionary where keys are attributes (strings),
                and values are dictionaries. The inner dictionaries have keys representing
                levels (strings) and values representing sets of proteome IDs (strings).
        """
        proteomes_by_level_by_attribute: Dict[str, Dict[str, Set[str]]] = {
            attribute: {} for attribute in self.attributes
        }
        for proteome_id in self.level_by_attribute_by_proteome_id:
            for attribute in self.attributes:
                level = self.level_by_attribute_by_proteome_id[proteome_id][attribute]
                if level not in proteomes_by_level_by_attribute[attribute]:
                    proteomes_by_level_by_attribute[attribute][level] = set()
                proteomes_by_level_by_attribute[attribute][level].add(proteome_id)
        return proteomes_by_level_by_attribute
