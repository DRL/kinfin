from collections import Counter
from typing import Any, Dict, List, Optional

from ete3 import Tree

from core.alo_collections import AloCollection
from core.logic import (
    add_taxid_attributes,
    parse_attributes_from_config_file,
    parse_attributes_from_json,
    parse_go_mapping,
    parse_ipr_mapping,
    parse_pfam_mapping,
    parse_tree_from_file,
)
from core.proteins import ProteinCollection
from core.utils import yield_file_lines


# cli
def parse_domains_from_functional_annotations_file(
    functional_annotation_f: str,
    proteinCollection: ProteinCollection,
    pfam_mapping: bool,
    ipr_mapping: bool,
    pfam_mapping_f: str,
    ipr_mapping_f: str,
    go_mapping_f: str,
) -> None:
    """
    Parse functional annotations from a file and populate ProteinCollection with parsed data.

    Parameters:
    - functional_annotation_f (str): Path to the functional annotation file.
    - proteinCollection (ProteinCollection): Instance of ProteinCollection class to store parsed data.
    - pfam_mapping (bool): Flag indicating whether to parse Pfam mappings.
    - ipr_mapping (bool): Flag indicating whether to parse InterPro mappings.
    - pfam_mapping_f (str): File path to the Pfam mapping file.
    - ipr_mapping_f (str): File path to the InterPro mapping file.
    - go_mapping_f (str): File path to the GO mapping file.

    Raises:
    - ValueError: If the functional annotation file lacks a header.

    Notes:
    - The function reads each line of the functional annotation file, parses relevant data,
      and populates the proteinCollection with domain annotations and GO terms.
    - It also optionally parses additional mappings (Pfam, InterPro, GO) based on provided flags.
    - Updates proteinCollection.functional_annotation_parsed and proteinCollection.domain_desc_by_id_by_source.
    """

    print(f"[STATUS] - Parsing {functional_annotation_f} ... this may take a while")

    for line in yield_file_lines(functional_annotation_f):
        temp: List[str] = line.split()
        if temp[0].startswith("#"):
            proteinCollection.domain_sources = temp[1:]

        else:
            if not proteinCollection.domain_sources:
                error_msg = f"[ERROR] - {functional_annotation_f} does not seem to have a header."  # fmt: skip
                raise ValueError(error_msg)

            domain_protein_id: str = temp.pop(0)
            go_terms: List[str] = []
            domain_counter_by_domain_source: Dict[str, Counter[str]] = {}
            for idx, field in enumerate(temp):
                if not field == "None":
                    domain_source: str = proteinCollection.domain_sources[idx]
                    domain_string: List[str] = field.split(";")
                    domain_counts_by_domain_id: Dict[str, int] = {}
                    for domain_id_count in domain_string:
                        domain_id: str
                        domain_count: int = 1
                        if domain_source == "GO":
                            domain_id = domain_id_count
                        else:
                            domain_id, domain_count_str = domain_id_count.rsplit(":", 2)  # fmt: skip
                            domain_count = int(domain_count_str)
                        domain_counts_by_domain_id[domain_id] = domain_count
                    domain_counter: Counter[str] = Counter(domain_counts_by_domain_id)  # fmt: skip
                    domain_counter_by_domain_source[domain_source] = domain_counter
            proteinCollection.add_annotation_to_protein(
                domain_protein_id=domain_protein_id,
                domain_counter_by_domain_source=domain_counter_by_domain_source,
                go_terms=go_terms,
            )

    proteinCollection.functional_annotation_parsed = True

    domain_desc_by_id_by_source: Dict[str, Dict[str, str]] = {}

    if pfam_mapping and "Pfam" in proteinCollection.domain_sources:
        domain_desc_by_id_by_source["Pfam"] = parse_pfam_mapping(pfam_mapping_f=pfam_mapping_f)  # fmt: skip

    if ipr_mapping and "IPR" in proteinCollection.domain_sources:
        domain_desc_by_id_by_source["IPR"] = parse_ipr_mapping(ipr_mapping_f=ipr_mapping_f)  # fmt: skip

    if go_mapping_f:
        domain_desc_by_id_by_source["GO"] = parse_go_mapping(go_mapping_f=go_mapping_f)  # fmt: skip

    proteinCollection.domain_desc_by_id_by_source = domain_desc_by_id_by_source


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


# api
def build_AloCollection_from_json(
    nodesdb_f: str,
    taxranks: List[str],
    json_list: List[Dict[str, str]],
    taxon_idx_mapping_file: str,
    tree_f: Optional[str],
):
    """
    Builds an AloCollection object from API input.

    Args:
        json_list List[Dict[str,str]]: JSON list of attributes.
        taxon_idx_mapping_file str: The path to the taxon-idx mapping file
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
