import json
import os
from collections import defaultdict
from typing import DefaultDict, Dict, List, Literal, Optional, Set, Tuple

import ete3
from ete3 import Tree, TreeNode

from core.utils import logger, progress, read_fasta_len, yield_file_lines


# common
def parse_nodesdb(filepath: str) -> Dict[str, Dict[str, str]]:
    """
    Parses the nodes database file.

    Args:
        filepath (str): The path to the nodes database file.

    Returns:
        Dict[str, Dict[str, str]]: A dictionary containing node information.
            Keys are node identifiers, and values are dictionaries with keys:

                - 'rank': The rank of the node.
                - 'name': The name of the node.
                - 'parent': The parent of the node.
    """
    logger.info(f"[STATUS] - Parsing nodesDB {filepath}")

    nodesdb: Dict[str, Dict[str, str]] = {}
    nodesdb_count = 0
    nodes_count = 0

    for line in yield_file_lines(filepath):
        if line.startswith("#"):
            nodesdb_count = int(line.lstrip("# nodes_count = ").rstrip("\n"))
        elif not line.strip():
            pass
        else:
            nodes_count += 1
            try:
                node, rank, name, parent = line.rstrip("\n").split("\t")
                nodesdb[node] = {"rank": rank, "name": name, "parent": parent}
            except Exception:
                pass
            if nodesdb_count:
                progress(nodes_count, 1000, nodesdb_count)
    return nodesdb


# cli
def get_lineage(
    taxid: str,
    nodesdb: Dict[str, Dict[str, str]],
    taxranks: List[str],
) -> Dict[str, str]:
    """
    Get the lineage of a taxonomic identifier.

    Args:
        taxid (str): The taxonomic identifier.
        nodesdb (Dict[str, Dict[str, str]]): A dictionary containing information about nodes.
        taxranks (List[str]): A list of taxonomic ranks to include in the lineage.

    Returns:
        Dict[str, str]: A dictionary containing the lineage information, with taxonomic ranks as keys
        and corresponding names as values.
    """
    lineage = {taxrank: "undef" for taxrank in taxranks}
    parent = ""
    node = taxid
    while parent != "1":
        taxrank = nodesdb[node]["rank"]
        name = nodesdb[node]["name"]
        parent = nodesdb[node]["parent"]
        if taxrank in taxranks:
            lineage[taxrank] = name
        node = parent
    return lineage


# cli
def parse_attributes_from_config_file(
    config_f: str,
) -> Tuple[Set[str], Dict[str, str], List[str], Dict[str, Dict[str, str]]]:
    """
    Parses attributes from a configuration file.

    Args:
        config_f (str): The path to the configuration file.

    Returns:
        Tuple[Set[str], Dict[str, str], List[str], Dict[str, Dict[str, str]]]: A tuple containing:
            - A set of proteome IDs.
            - A dictionary mapping species IDs to proteome IDs.
            - A list of attributes.
            - A dictionary mapping proteome IDs to dictionaries, where each inner dictionary
              maps attributes to their corresponding levels.

    Raises:
        FileNotFoundError: If the specified configuration file is not found.
        ValueError: If there are errors in the configuration file format or content.

    Note:
        - The configuration file is expected to have a header line starting with '#',
          where the first element is 'IDX' and the second element is 'TAXON'.
        - Each subsequent non-empty line in the configuration file should contain
          comma-separated values corresponding to the attributes defined in the header line.
        - The 'TAXON' attribute is expected to be unique for each line.
    """

    logger.info(f"[STATUS] - Parsing config file: {config_f} ...")
    attributes: List[str] = []
    level_by_attribute_by_proteome_id: Dict[str, Dict[str, str]] = {}
    proteomes: Set[str] = set()
    proteome_id_by_species_id: Dict[str, str] = {}

    for line in yield_file_lines(config_f):
        if line.startswith("#"):
            if not attributes:
                attributes = [x.strip() for x in line.lstrip("#").split(",")]
                if not attributes[0] == "IDX" or not attributes[1] == "TAXON":
                    error_msg = f"[ERROR] - First/second element have to be IDX/TAXON.\n\t{attributes}"
                    raise ValueError(error_msg)
            else:
                # accounts for SpeciesIDs that are commented out for Orthofinder
                pass

        elif line.strip():
            temp = line.split(",")

            if not len(temp) == len(attributes):
                error_msg = f"[ERROR] - number of columns in line differs from header\n\t{attributes}\n\t{temp}"
                raise ValueError(error_msg)

            if temp[1] in proteomes:
                error_msg = f"[ERROR] - 'TAXON' should be unique. {temp[0]} was encountered multiple times"
                raise ValueError(error_msg)

            species_id = temp[0]
            proteome_id = temp[1]
            proteomes.add(proteome_id)
            proteome_id_by_species_id[species_id] = proteome_id

            level_by_attribute_by_proteome_id[proteome_id] = {
                attribute: level for attribute, level in zip(attributes, temp)
            }
            level_by_attribute_by_proteome_id[proteome_id]["all"] = "all"
        else:
            pass
    attributes.insert(0, "all")  # append to front
    return (
        proteomes,
        proteome_id_by_species_id,
        attributes,
        level_by_attribute_by_proteome_id,
    )


# common
def add_taxid_attributes(
    nodesdb_f: str,
    taxranks: List[str],
    attributes: List[str],
    level_by_attribute_by_proteome_id: Dict[str, Dict[str, str]],
) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """
    Adds taxonomic attributes to the dictionary of attributes indexed by proteome ID.

    Parameters:

        - nodesdb_f (str): File path to the nodes database.
        - taxranks (List[str]): List of taxonomic ranks to be included as attributes.
        - attributes (List[str]): List of existing attributes.
        - level_by_attribute_by_proteome_id (Dict[str, Dict[str, str]]): Dictionary where keys
            are proteome IDs and values are dictionaries of attributes for each proteome ID,
            including at least the "TAXID" attribute.

    Returns:
        Tuple[List[str], Dict[str, Dict[str, str]]]: A tuple containing:

            - Updated list of attributes with taxonomic ranks added and "TAXID" removed.
            - Updated dictionary of attributes indexed by proteome ID, with taxonomic attributes added and "TAXID" removed.
    """
    NODESDB = parse_nodesdb(str(nodesdb_f))
    for proteome_id in level_by_attribute_by_proteome_id:
        taxid = level_by_attribute_by_proteome_id[proteome_id]["TAXID"]
        lineage = get_lineage(taxid=taxid, nodesdb=NODESDB, taxranks=taxranks)

        # add lineage attribute/levels
        for taxrank in taxranks:
            level_by_attribute_by_proteome_id[proteome_id][taxrank] = lineage[taxrank]

        # remove taxid-levels
        del level_by_attribute_by_proteome_id[proteome_id]["TAXID"]

    # remove taxid-attribute
    attributes.remove("TAXID")

    # add taxranks to rank
    for taxrank in taxranks:
        attributes.append(taxrank)

    return attributes, level_by_attribute_by_proteome_id


# cli
def parse_tree_from_file(
    tree_f: str,
    outgroups: List[str],
) -> tuple[Tree, Dict[frozenset[str], str]]:
    """
    Parse a phylogenetic tree from nwk file and set specified outgroups.

    Args:
        tree_f (str): Path to the nwk tree file.
        outgroups (List[str]): List of outgroup taxa names.

    Returns:
        tuple[ete3.Tree, Dict[str, int]]: A tuple containing the parsed phylogenetic tree
            and a dictionary mapping proteome IDs to node indices.
    """
    logger.info(f"[STATUS] - Parsing Tree file : {tree_f} ...")
    tree_ete: TreeNode = ete3.Tree(tree_f)
    if len(outgroups) > 1:
        outgroup_node: TreeNode = tree_ete.get_common_ancestor(outgroups)  # type: ignore
        try:
            logger.info(f"[STATUS] - Setting LCA of {", ".join(outgroups)} as outgroup : ...")
            tree_ete.set_outgroup(outgroup_node)  # type: ignore
        except ete3.coretype.tree.TreeError:  # type: ignore
            logger.info("[STATUS] - Tree seems to be rooted already : ...")
    else:
        logger.info(f"[STATUS] - Setting {",".join(outgroups)} as outgroup : ...")
        tree_ete.set_outgroup(outgroups[0])  # type: ignore
    logger.info(tree_ete)
    node_idx_by_proteome_ids: Dict[frozenset[str], str] = {}
    for idx, node in enumerate(tree_ete.traverse("levelorder")):  # type: ignore
        proteome_ids = frozenset([leaf.name for leaf in node])
        if not node.name:
            node.add_features(
                name=f"n{idx}",
                nodetype="node",
                proteome_ids=proteome_ids,
                apomorphic_cluster_counts={"singletons": 0, "non_singletons": 0},
                synapomorphic_cluster_counts={
                    "complete_presence": 0,
                    "partial_absence": 0,
                },
                synapomorphic_cluster_strings=[],
                counts={"specific": 0, "shared": 0, "absent": 0, "singleton": 0},
            )
        else:
            node.add_features(
                nodetype="tip",
                proteome_ids=proteome_ids,
                apomorphic_cluster_counts={"singletons": 0, "non_singletons": 0},
                synapomorphic_cluster_counts={
                    "complete_presence": 0,
                    "partial_absence": 0,
                },
                synapomorphic_cluster_strings=[],
                counts={"specific": 0, "shared": 0, "absent": 0, "singleton": 0},
            )
        node_idx_by_proteome_ids[proteome_ids] = node.name
    return tree_ete, node_idx_by_proteome_ids


# api
def parse_attributes_from_json(
    json_list: List[Dict[str, str]],
    taxon_idx_mapping_file: str,
) -> Tuple[Set[str], Dict[str, str], List[str], Dict[str, Dict[str, str]]]:
    """
    Parses attributes from a JSON list.

    Args:
        json_list List[Dict[str,str]]: JSON list of attributes.
        taxon_idx_mapping_file str: The path to the taxon-idx mapping file

    Returns:
        Tuple[Set[str], Dict[str, str], List[str], Dict[str, Dict[str, str]]]: A tuple containing:
            - A set of proteome IDs.
            - A dictionary mapping species IDs to proteome IDs.
            - A list of attributes.
            - A dictionary mapping proteome IDs to dictionaries, where each inner dictionary
              maps attributes to their corresponding levels.

    Raises:
        FileNotFoundError: If the specified configuration file is not found.
        ValueError: If there are errors in the configuration file format or content.

    Note:
        - The configuration file is expected to have a header line starting with '#',
          where the first element is 'IDX' and the second element is 'TAXON'.
        - Each subsequent non-empty line in the configuration file should contain
          comma-separated values corresponding to the attributes defined in the header line.
        - The 'TAXON' attribute is expected to be unique for each line.
    """

    logger.info("[STATUS] - Parsing JSON list...")
    attributes: List[str] = []
    level_by_attribute_by_proteome_id: Dict[str, Dict[str, str]] = {}
    proteomes: Set[str] = set()
    proteome_id_by_species_id: Dict[str, str] = {}

    attributes = list(json_list[0].keys())
    attributes.insert(0, "IDX")

    with open(taxon_idx_mapping_file, "r") as f:
        taxon_idx_mapping = json.load(f)

    attributes.insert(0, "all")

    for entry in json_list:
        proteome_id = entry["TAXON"]
        species_id = taxon_idx_mapping[proteome_id]
        proteomes.add(proteome_id)
        proteome_id_by_species_id[species_id] = proteome_id

        level_by_attribute_by_proteome_id[proteome_id] = {
            attribute: entry.get(attribute, "") for attribute in attributes[1:]
        }
        level_by_attribute_by_proteome_id[proteome_id]["IDX"] = proteome_id
        level_by_attribute_by_proteome_id[proteome_id]["all"] = "all"
    attributes.insert(0, "all")
    return (
        proteomes,
        proteome_id_by_species_id,
        attributes,
        level_by_attribute_by_proteome_id,
    )


def parse_fasta_dir(species_ids_f: str, fasta_dir: str) -> Dict[str, int]:
    """
    Parse a species IDs file to retrieve fasta file names and then calculate
    lengths of sequences from corresponding FASTA files.

    Args:
    - species_ids_f (str): Path to the species IDs file, where each line contains
      an index and a corresponding FASTA file name separated by ': '.
    - fasta_dir (str): Directory path where the FASTA files are located.

    Returns:
    - Dict[str, int]: A dictionary mapping header strings (protein IDs) to their
      corresponding sequence lengths extracted from the FASTA files.
    """
    logger.info("[STATUS] - Parsing FASTAs ...")
    fasta_file_by_species_id: Dict[str, str] = {}

    for line in yield_file_lines(species_ids_f):
        if not line.startswith("#"):
            idx, fasta = line.split(": ")
            fasta_file_by_species_id[idx] = fasta

    fasta_len_by_protein_id: Dict[str, int] = {}
    for _, fasta_f in list(fasta_file_by_species_id.items()):
        fasta_f = os.path.join(fasta_dir, fasta_f)

        for header, length in read_fasta_len(fasta_f):
            fasta_len_by_protein_id[header] = length

    return fasta_len_by_protein_id


def parse_pfam_mapping(pfam_mapping_f: str) -> Dict[str, str]:
    """
    Parse a PFAM mapping file to create a dictionary mapping PFAM domain IDs to their descriptions.

    Args:
    - pfam_mapping_f (str): Path to the PFAM mapping file, where each line contains tab-separated values
      with the domain ID in the first column and its description in the fifth column.

    Returns:
    - Dict[str, str]: A dictionary mapping PFAM domain IDs to their corresponding descriptions.

    Raises:
    - ValueError: If conflicting descriptions are found for the same domain ID.
    """
    logger.info(f"[STATUS] - Parsing {pfam_mapping_f} ... ")

    pfam_mapping_dict: Dict[str, str] = {}
    for line in yield_file_lines(pfam_mapping_f):
        temp: List[str] = line.split("\t")
        domain_id: str = temp[0]
        domain_desc: str = temp[4]
        if domain_id not in pfam_mapping_dict:
            pfam_mapping_dict[domain_id] = domain_desc
        else:
            if not domain_desc == pfam_mapping_dict[domain_id]:
                error_msg = f"[ERROR] : Conflicting descriptions for {domain_id}"
                raise ValueError(error_msg)

    return pfam_mapping_dict


def parse_ipr_mapping(ipr_mapping_f: str) -> Dict[str, str]:
    """
    Parse an InterPro (IPR) mapping file to create a dictionary mapping InterPro IDs to their descriptions.

    Args:
    - ipr_mapping_f (str): Path to the InterPro mapping file, where each line contains an InterPro ID and its description.
      Lines starting with "Active_site" are skipped as they are not relevant to mapping.

    Returns:
    - Dict[str, str]: A dictionary mapping InterPro IDs to their corresponding descriptions.

    Raises:
    - ValueError: If conflicting descriptions are found for the same InterPro ID.
    """
    logger.info(f"[STATUS] - Parsing {ipr_mapping_f} ... ")

    ipr_mapping_dict: Dict[str, str] = {}
    for line in yield_file_lines(ipr_mapping_f):
        if not line.startswith("Active_site"):
            temp: List[str] = line.split()
            ipr_id: str = temp[0]
            ipr_desc: str = " ".join(temp[1:])
            if ipr_id not in ipr_mapping_dict:
                ipr_mapping_dict[ipr_id] = ipr_desc
            else:
                if not ipr_desc == ipr_mapping_dict[ipr_id]:
                    error_msg = f"[ERROR] : Conflicting descriptions for {ipr_id}"
                    raise ValueError(error_msg)
    return ipr_mapping_dict


def parse_go_mapping(go_mapping_f: str) -> Dict[str, str]:
    """
    Parse a Gene Ontology (GO) mapping file to create a dictionary mapping GO IDs to their descriptions.

    Args:
    - go_mapping_f (str): Path to the GO mapping file, where each line contains a GO ID and its description.
      Lines starting with '!' are skipped as they are comments.

    Returns:
    - Dict[str, str]: A dictionary mapping GO IDs (without 'GO:' prefix) to their corresponding descriptions.

    Raises:
    - ValueError: If conflicting descriptions are found for the same GO ID.
    """
    logger.info(f"[STATUS] - Parsing {go_mapping_f} ... ")
    go_mapping_dict: Dict[str, str] = {}
    for line in yield_file_lines(go_mapping_f):
        if not line.startswith("!"):
            temp: List[str] = line.replace(" > ", "|").split("|")
            go_string: List[str] = temp[1].split(";")
            go_desc, go_id = go_string[0].replace("GO:", ""), go_string[1].lstrip(" ")

            if go_id not in go_mapping_dict:
                go_mapping_dict[go_id] = go_desc
            else:
                if not go_desc == go_mapping_dict[go_id]:
                    error_msg = f"[ERROR] : Conflicting descriptions for {go_id}"
                    raise ValueError(error_msg)
    return go_mapping_dict


def compute_protein_ids_by_proteome(
    proteomes_by_protein_id: Dict[str, str]
) -> DefaultDict[str, Set[str]]:
    """
    Compute protein IDs grouped by proteome IDs.

    Args:
        proteomes_by_protein_id (Dict[str, str]): A dictionary mapping protein IDs to proteome IDs.

    Returns:
        DefaultDict[str, Set[str]]: A defaultdict where keys are proteome IDs and values are sets
        of protein IDs belonging to each proteome ID.
    """
    protein_ids_by_proteome_id: DefaultDict[str, Set[str]] = defaultdict(set)
    for protein_id, proteome_id in list(proteomes_by_protein_id.items()):
        protein_ids_by_proteome_id[proteome_id].add(protein_id)
    return protein_ids_by_proteome_id


# common
def get_attribute_cluster_type(
    singleton,
    implicit_protein_ids_by_proteome_id_by_level,
) -> Literal["singleton", "shared", "specific"]:
    """
    Determines the type of cluster based on the parameters.

    Parameters:
    - singleton: A boolean indicating whether the cluster is a singleton.
    - implicit_protein_ids_by_proteome_id_by_level: A dictionary representing protein ids 
      grouped by proteome id at different levels.

    Returns:
    - One of the following strings:
      - "singleton": If `singleton` is True.
      - "shared": If there are protein ids grouped under multiple proteome ids.
      - "specific": If there is only one proteome id with protein ids.

    """
    if singleton:
        return "singleton"
    else:
        if len(implicit_protein_ids_by_proteome_id_by_level) > 1:
            return "shared"
        else:
            return "specific"


def get_ALO_cluster_cardinality(
    ALO_proteome_counts_in_cluster: List[int],
    fuzzy_range: Set[int],
    fuzzy_count: int = 1,
    fuzzy_fraction: float = 0.75,
) -> Optional[str]:
    """
    Determine the cardinality type of a cluster based on ALO proteome counts.

    Args:
        ALO_proteome_counts_in_cluster (List[int]): List of ALO proteome counts in the cluster.
        fuzzy_range (Set[int]): Set of integers representing the range of fuzzy counts.
        fuzzy_count (int, optional): Specific count considered as fuzzy. Default is 1.
        fuzzy_fraction (float, optional): Fraction threshold for considering a cluster as 'fuzzy'. Default is 0.75.

    Returns:
        Optional[str]: Returns "true" (str) if all counts are 1, "fuzzy" (str) if the cluster meets fuzzy criteria,
                       and None otherwise.
    """
    if len(ALO_proteome_counts_in_cluster) > 2:
        length = len(ALO_proteome_counts_in_cluster)
        if all(count == 1 for count in ALO_proteome_counts_in_cluster):
            return "true"
        else:
            fuzzycount_count = len(
                [
                    ALO_proteome_counts
                    for ALO_proteome_counts in ALO_proteome_counts_in_cluster
                    if ALO_proteome_counts == fuzzy_count
                ]
            )

            fuzzyrange_count = len(
                [
                    ALO_proteome_counts
                    for ALO_proteome_counts in ALO_proteome_counts_in_cluster
                    if ALO_proteome_counts in fuzzy_range
                ]
            )

            fuzzy_fr = fuzzycount_count / length

            if fuzzy_fr >= fuzzy_fraction:
                if fuzzycount_count + fuzzyrange_count == length:
                    return "fuzzy"

    return None

