from typing import Dict, List, Set, Tuple

from core.utils import progress, yield_file_lines


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
    print(f"[STATUS] - Parsing nodesDB {filepath}")

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

    print(f"[STATUS] - Parsing config file: {config_f} ...")
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
