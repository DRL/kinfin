from typing import Dict, List, Optional, Set, Tuple


class ServeArgs:
    def __init__(self, port: int = 8000):
        self.port = port


class InputData:
    def __init__(
        self,
        nodesdb_f: str,
        pfam_mapping_f: str,
        ipr_mapping_f: str,
        go_mapping_f: str,
        cluster_file: str,
        config_data: List[Dict[str, str]] | str,
        sequence_ids_file: str,
        species_ids_file: Optional[str] = None,
        functional_annotation_f: Optional[str] = None,
        fasta_dir: Optional[str] = None,
        tree_file: Optional[str] = None,
        output_path: Optional[str] = None,
        infer_singletons: Optional[bool] = False,
        plot_tree: bool = False,
        min_proteomes: int = 2,
        test: str = "mannwhitneyu",
        taxranks: List[str] = ["phylum", "order", "genus"],
        repetitions: int = 30,
        fuzzy_count: int = 1,
        fuzzy_fraction: float = 0.75,
        fuzzy_range: Set[int] = set([x for x in range(0, 20 + 1) if not x == 1]),
        fontsize: int = 18,
        plotsize: Tuple[float, float] = (24, 12),
        plot_format: str = "pdf",
        taxon_idx_mapping_file: Optional[str] = None,
    ):
        self.cluster_f = cluster_file
        self.config_data = config_data
        self.sequence_ids_f = sequence_ids_file
        self.species_ids_f = species_ids_file
        self.tree_f = tree_file
        self.functional_annotation_f = functional_annotation_f
        self.taxon_idx_mapping_file = taxon_idx_mapping_file

        self.nodesdb_f = nodesdb_f
        self.pfam_mapping_f = pfam_mapping_f
        self.ipr_mapping_f = ipr_mapping_f
        self.go_mapping_f = go_mapping_f

        self.test = test
        self.plot_tree = plot_tree
        self.fasta_dir = fasta_dir
        self.output_path = output_path
        self.infer_singletons = infer_singletons
        self.fuzzy_count = fuzzy_count
        self.fuzzy_fraction = fuzzy_fraction
        self.fuzzy_range = fuzzy_range
        self.repetitions = repetitions
        self.min_proteomes = min_proteomes
        self.plot_format = plot_format
        self.fontsize = fontsize
        self.taxranks = taxranks
        self.plotsize = plotsize

        self.pfam_mapping = True
        self.ipr_mapping = True
