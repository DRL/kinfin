from typing import Dict, List, Set, Tuple

from core.utils import yield_file_lines


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
