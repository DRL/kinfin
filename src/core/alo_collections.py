import os
from typing import Any, Dict, List, Optional, Set, Tuple

import ete3
import matplotlib as mat
import numpy as np
from ete3 import Tree

from core.alo import AttributeLevel
from core.config import ATTRIBUTE_RESERVED

mat.use("agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, NullFormatter

plt.style.use("ggplot")
mat.rc("ytick", labelsize=20)
mat.rc("xtick", labelsize=20)
axis_font = {"size": "20"}
mat.rcParams.update({"font.size": 22})


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
        self.fastas_parsed: bool = False
        self.ALO_by_level_by_attribute = self.create_ALOs()

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

    def create_ALOs(self) -> Dict[str, Dict[str, Optional[AttributeLevel]]]:
        """
        Creates Attribute Level Objects (ALOs) for each attribute and level based on
        proteome IDs.

        Returns:
            Dict[str, Dict[str, Optional[AttributeLevel]]]:
                A dictionary where each key is an attribute name (str),
                and the corresponding value is a dictionary mapping level (str)
                to an AttributeLevel instance or None.
        """
        ALO_by_level_by_attribute: Dict[str, Dict[str, Optional[AttributeLevel]]] = { attribute: {} for attribute in self.attributes }  # fmt:skip
        for attribute in self.proteome_ids_by_level_by_attribute:
            for level in self.proteome_ids_by_level_by_attribute[attribute]:
                proteome_ids = self.proteome_ids_by_level_by_attribute[attribute][level]
                ALO = AttributeLevel(
                    #
                    attribute=attribute,
                    level=level,
                    proteomes=proteome_ids,
                )
                if level not in ALO_by_level_by_attribute[attribute]:
                    ALO_by_level_by_attribute[attribute][level] = None
                ALO_by_level_by_attribute[attribute][level] = ALO
        return ALO_by_level_by_attribute

    def generate_header_for_node(self, node: ete3.TreeNode, dirs: Dict[str, str]):
        """
        Generates a header image for a given node of a tree with specified statistics.

        Args:
            node (ete3.TreeNode): The TreeNode object representing the node for which the header is generated.
            dirs (Dict[str, str]): A dictionary containing directory paths, including 'tree_headers' where the header image will be saved.

        Returns:
            str: File path to the generated header image.

        Notes:
            - The method generates a header image in PNG format displaying various statistics (apomorphies and synapomorphies) for the given tree node.
            - The statistics include counts of singletons, non-singletons, complete presence synapomorphies, and partial absence synapomorphies.
            - The generated image is saved in the specified directory under 'tree_headers' with the node's name as the filename.

        Raises:
            Any exceptions that might occur during file saving or table rendering.
        """

        node_header_f = os.path.join(dirs["tree_headers"], f"{node.name}.header.png")
        data = []
        data.append(
            (
                "Apomorphies (size=1)",
                "{:,}".format(
                    node.apomorphic_cluster_counts["singletons"]  # type:ignore
                ),
            )
        )
        data.append(
            (
                "Apomorphies (size>1)",
                "{:,}".format(
                    node.apomorphic_cluster_counts["non_singletons"]  # type:ignore
                ),
            )
        )
        data.append(
            (
                "Synapomorphies (all)",
                "{:,}".format(
                    node.synapomorphic_cluster_counts[  # type:ignore
                        "complete_presence"
                    ]
                    + node.synapomorphic_cluster_counts[  # type:ignore
                        "partial_absence"
                    ]
                ),
            )
        )
        data.append(
            (
                "Synapomorphies (cov=100%)",
                "{:,}".format(
                    node.synapomorphic_cluster_counts[  # type:ignore
                        "complete_presence"
                    ]
                ),
            )
        )
        data.append(
            (
                "Synapomorphies (cov<100%)",
                "{:,}".format(
                    node.synapomorphic_cluster_counts["partial_absence"]  # type:ignore
                ),
            )
        )
        col_labels = ("Type", "Count")
        fig, ax = plt.subplots(figsize=(2, 0.5))
        ax.set_facecolor("white")
        table = ax.table(
            cellText=data,
            colLabels=col_labels,
            loc="bottom",
            fontsize=24,
            colLoc="center",
            rowLoc="right",
            edges="",
        )
        table.set_fontsize(24)
        table.scale(2, 1)
        for key, cell in list(table.get_celld().items()):
            row, col = key
            cell._text.set_color("grey")
            if row > 0:
                cell.set_edgecolor("darkgrey")
                cell.visible_edges = "T"
            else:
                cell.set_edgecolor("darkgrey")
                cell.visible_edges = "B"
            if row == len(data) - 2:
                cell.set_edgecolor("darkgrey")
                cell.visible_edges = "T"
        ax.axis("tight")
        ax.axis("off")
        print(f"[STATUS]\t- Plotting {node_header_f}")
        fig.savefig(node_header_f, pad=0, bbox_inches="tight", format="png")
        plt.close()
        return node_header_f

    def generate_chart_for_node(
        self,
        node,
        dirs: Dict[str, str],
        plot_format: str,
        fontsize: int,
    ) -> Optional[str]:
        """
        Generate and save a histogram chart for a given node's synapomorphies.

        Args:
        - node: The node object containing synapomorphic cluster strings.
        - dirs: A dictionary containing directory paths, specifically 'tree_charts' for saving charts.
        - plot_format: The format in which to save the chart ('png' or 'pdf').
        - fontsize: Font size for axis labels and ticks.

        Returns:
        - Optional[str]: Path to the saved chart file if successful, None otherwise.
        """

        proteome_coverages = []
        for synapomorphic_cluster_string in node.synapomorphic_cluster_strings:
            proteome_coverages.append(float(synapomorphic_cluster_string[3]))
        if proteome_coverages:
            chart_f = os.path.join(dirs["tree_charts"], f"{node.name}.barchart.png")
            f, ax = plt.subplots(figsize=(3.0, 3.0))
            ax.set_facecolor("white")
            x_values = np.array(proteome_coverages)
            ax.hist(
                x_values,
                histtype="stepfilled",
                align="mid",
                bins=np.arange(0.0, 1.0 + 0.1, 0.1),
            )
            ax.set_xlim(-0.1, 1.1)
            for tick in ax.xaxis.get_majorticklabels():
                tick.set_fontsize(fontsize - 2)
                tick.set_rotation("vertical")
            for tick in ax.yaxis.get_majorticklabels():
                tick.set_fontsize(fontsize - 2)
            ax.set_frame_on(False)
            ax.xaxis.grid(True, linewidth=1, which="major", color="lightgrey")
            ax.yaxis.grid(True, linewidth=1, which="major", color="lightgrey")
            f.suptitle("Synapomorphies", y=1.1)
            ax.set_ylabel("Count", fontsize=fontsize)
            ax.set_xlabel("Proteome coverage", fontsize=fontsize)
            print(f"[STATUS]\t- Plotting {chart_f}")
            f.savefig(chart_f, bbox_inches="tight", format="png")
            if plot_format == "pdf":
                pdf_chart_f = os.path.join(
                    dirs["tree_charts"], f"{node.name}.barchart.pdf"
                )
                print(f"[STATUS]\t- Plotting {pdf_chart_f}")
                f.savefig(pdf_chart_f, bbox_inches="tight", format="pdf")
            plt.close()
            return chart_f

    def plot_text_tree(self, dirs: Dict[str, str]) -> None:
        """
        Plot and save the textual representation of the tree.

        This method uses the `tree_ete` attribute of the class to generate and save
        both a Newick format (.nwk) and a text format (.txt) representation of the tree.

        Args:
        - dirs: A dictionary containing directory paths, specifically 'tree' for saving tree files.

        Returns:
        - None
        """
        if self.tree_ete:
            tree_nwk_f = os.path.join(dirs["tree"], "tree.nwk")
            self.tree_ete.write(format=1, outfile=tree_nwk_f)
            tree_txt_f = os.path.join(dirs["tree"], "tree.txt")
            with open(tree_txt_f, "w") as tree_txt_fh:
                tree_txt_fh.write(f"{self.tree_ete.get_ascii(show_internal=True, compact=False)}\n")  # fmt:skip
