from typing import Any, Dict, List, Optional

from ete3 import Tree

from core.alo_collections import AloCollection
from core.logic import (
    add_taxid_attributes,
    parse_attributes_from_config_file,
    parse_attributes_from_json,
    parse_tree_from_file,
)


# cli
def build_AloCollection(
    config_f: str,
    nodesdb_f: str,
    taxranks: List[str],
    tree_f: Optional[str],
) -> AloCollection:
    """
    Builds an AloCollection object from command-line interface (CLI) inputs.

    Args:
        config_f (str): Path to the configuration file containing proteome attributes.
        nodesdb_f (str): Path to the nodes database file for inferring taxonomic ranks.
        taxranks (List[str]): List of taxonomic ranks to be inferred.
        tree_f (Optional[str]): Path to the tree file. If provided, ALOs are added from the tree.

    Returns:
        AloCollection: An instance of the AloCollection class containing parsed data.
    """
    (
        proteomes,
        proteome_id_by_species_id,
        attributes,
        level_by_attribute_by_proteome_id,
    ) = parse_attributes_from_config_file(config_f)

    # Add taxonomy if needed
    if "TAXID" in set(attributes):
        print("[STATUS] - Attribute 'TAXID' found, inferring taxonomic ranks from nodesDB")  # fmt: skip
        attributes, level_by_attribute_by_proteome_id = add_taxid_attributes(
            attributes=attributes,
            level_by_attribute_by_proteome_id=level_by_attribute_by_proteome_id,
            nodesdb_f=nodesdb_f,
            taxranks=taxranks,
        )

    # Add ALOs from tree if provided
    tree_ete: Optional[Tree] = None
    node_idx_by_proteome_ids: Optional[Dict[Any, Any]] = None
    if tree_f:
        outgroups: List[str] = []
        if "OUT" not in attributes:
            error_msg = "[ERROR] - Please specify one of more outgroup taxa"  # fmt: skip
            ValueError(error_msg)
        outgroups = [
            proteome_id
            for proteome_id in proteomes
            if level_by_attribute_by_proteome_id[proteome_id]["OUT"] == "1"
        ]
        tree_ete, node_idx_by_proteome_ids = parse_tree_from_file(tree_f, outgroups)

    print("[STATUS] - Building AloCollection ...")
    return AloCollection(
        proteomes=proteomes,
        attributes=attributes,
        proteome_id_by_species_id=proteome_id_by_species_id,
        level_by_attribute_by_proteome_id=level_by_attribute_by_proteome_id,
        node_idx_by_proteome_ids=node_idx_by_proteome_ids,
        tree_ete=tree_ete,
    )


def build_ALoCollection_from_json(
    nodesdb_f: str,
    taxranks: List[str],
    json_list: List[Dict[str, str]],
    taxon_idx_mapping_file: str,
    tree_f: Optional[str],
):
    (
        proteomes,
        proteome_id_by_species_id,
        attributes,
        level_by_attribute_by_proteome_id,
    ) = parse_attributes_from_json(
        json_list=json_list,
        taxon_idx_mapping_file=taxon_idx_mapping_file,
    )

    # Add taxonomy if needed
    if "TAXID" in set(attributes):
        print("[STATUS] - Attribute 'TAXID' found, inferring taxonomic ranks from nodesDB")  # fmt: skip
        attributes, level_by_attribute_by_proteome_id = add_taxid_attributes(
            attributes=attributes,
            level_by_attribute_by_proteome_id=level_by_attribute_by_proteome_id,
            nodesdb_f=nodesdb_f,
            taxranks=taxranks,
        )
    # Add ALOs from tree if provided
    tree_ete: Optional[Tree] = None
    node_idx_by_proteome_ids: Optional[Dict[Any, Any]] = None
    if tree_f:
        outgroups: List[str] = []
        if "OUT" not in attributes:
            error_msg = "[ERROR] - Please specify one of more outgroup taxa"  # fmt: skip
            ValueError(error_msg)
        outgroups = [
            proteome_id
            for proteome_id in proteomes
            if level_by_attribute_by_proteome_id[proteome_id]["OUT"] == "1"
        ]
        tree_ete, node_idx_by_proteome_ids = parse_tree_from_file(tree_f, outgroups)

    print("[STATUS] - Building AloCollection ...")
    return AloCollection(
        proteomes=proteomes,
        attributes=attributes,
        proteome_id_by_species_id=proteome_id_by_species_id,
        level_by_attribute_by_proteome_id=level_by_attribute_by_proteome_id,
        node_idx_by_proteome_ids=node_idx_by_proteome_ids,
        tree_ete=tree_ete,
    )
