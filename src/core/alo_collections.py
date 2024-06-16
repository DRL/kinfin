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

    def plot_tree(
        self,
        header_f_by_node_name,
        charts_f_by_node_name,
        dirs: Dict[str, str],
    ) -> None:
        """
        Plot and save a tree visualization with custom header and chart images for nodes.

        This method uses the `self.tree_ete` attribute of the class to visualize the tree
        in a hierarchical manner, with customized header and chart images for each node.

        Args:
        - header_f_by_node_name: Dictionary mapping node names to header image file paths (must be PNG).
        - charts_f_by_node_name: Dictionary mapping node names to chart image file paths (must be PNG).
        - dirs: A dictionary containing directory paths, specifically 'tree' for saving the tree PDF.

        Returns:
        - None
        """
        tree_f = os.path.join(
            dirs["tree"], "tree.pdf"
        )  # must be PDF! (otherwise it breaks)
        style = ete3.NodeStyle()
        style["vt_line_width"] = 5
        style["hz_line_width"] = 5
        style["fgcolor"] = "darkgrey"
        for node in self.tree_ete.traverse("levelorder"):  # type: ignore
            node.set_style(style)
            if header_f_by_node_name[node.name]:
                node_header_face = ete3.faces.ImgFace(
                    header_f_by_node_name[node.name]
                )  # must be PNG! (ETE can't do PDF Faces)
                node.add_face(node_header_face, column=0, position="branch-top")
            if charts_f_by_node_name[node.name]:
                node_chart_face = ete3.faces.ImgFace(
                    charts_f_by_node_name[node.name]
                )  # must be PNG! (ETE can't do PDF Faces)
                node.add_face(node_chart_face, column=0, position="branch-bottom")
            node_name_face = ete3.TextFace(node.name, fsize=64)
            node.img_style["size"] = 10
            node.img_style["shape"] = "sphere"
            node.img_style["fgcolor"] = "black"
            if not node.is_leaf():
                node.add_face(node_name_face, column=0, position="branch-right")
            node.add_face(node_name_face, column=0, position="aligned")
        ts = ete3.TreeStyle()
        ts.draw_guiding_lines = True
        ts.show_scale = False
        ts.show_leaf_name = False
        ts.allow_face_overlap = True
        ts.guiding_lines_color = "lightgrey"
        print(f"[STATUS] - Writing tree {tree_f}... ")
        self.tree_ete.render(tree_f, dpi=600, h=1189, units="mm", tree_style=ts)  # type: ignore

    def write_tree(
        self,
        dirs: Dict[str, str],
        render_tree: bool,
        plot_format: str,
        fontsize: int,
    ) -> None:
        """
        Write tree data to files and optionally render a graphical tree representation.

        This method generates and saves various metrics and data related to the tree structure,
        including node statistics and cluster metrics. It can also render a graphical tree
        representation if specified.

        Args:
        - dirs: A dictionary containing directory paths, including 'tree' for saving tree-related files.
        - render_tree: Boolean flag indicating whether to render a graphical tree representation.
        - plot_format: Format for saving plots ('png', 'pdf', etc.).
        - fontsize: Font size used for plotting.

        Returns:
        - None
        """
        if self.tree_ete:
            print("[STATUS] - Writing data for tree ... ")
            # Node stats
            node_stats_f = os.path.join(dirs["tree"], "tree.node_metrics.txt")
            node_stats_header: List[str] = []
            node_stats_header.append("#nodeID")
            node_stats_header.append("taxon_specific_apomorphies_singletons")
            node_stats_header.append("taxon_specific_apomorphies_non_singletons")
            node_stats_header.append("node_specific_synapomorphies_total")
            node_stats_header.append("node_specific_synapomorphies_complete_presence")
            node_stats_header.append("node_specific_synapomorphies_partial_absence")
            node_stats_header.append("proteome_count")
            node_stats: List[str] = []
            node_stats.append("\t".join(node_stats_header))
            # Cluster node stats
            node_clusters_f = os.path.join(dirs["tree"], "tree.cluster_metrics.txt")
            node_clusters_header = []
            node_clusters_header.append("#clusterID")
            node_clusters_header.append("nodeID")
            node_clusters_header.append("synapomorphy_type")
            node_clusters_header.append("node_taxon_coverage")
            node_clusters_header.append("children_coverage")
            node_clusters_header.append("node_taxa_present")
            node_clusters = []
            node_clusters.append("\t".join(node_clusters_header))
            # header_f_by_node_name
            header_f_by_node_name = {}
            charts_f_by_node_name = {}
            for node in self.tree_ete.traverse("levelorder"):  # type: ignore
                for synapomorphic_cluster_string in node.synapomorphic_cluster_strings:  # type: ignore
                    node_clusters.append(
                        "\t".join(
                            [
                                str(string)
                                for string in list(synapomorphic_cluster_string)
                            ]
                        )
                    )
                node_stats_line = []
                node_stats_line.append(node.name)
                node_stats_line.append(node.apomorphic_cluster_counts["singletons"])  # type: ignore
                node_stats_line.append(node.apomorphic_cluster_counts["non_singletons"])  # type: ignore
                node_stats_line.append(node.synapomorphic_cluster_counts["complete_presence"] + node.synapomorphic_cluster_counts["partial_absence"])  # type: ignore
                node_stats_line.append(node.synapomorphic_cluster_counts["complete_presence"])  # type: ignore
                node_stats_line.append(node.synapomorphic_cluster_counts["partial_absence"])  # type: ignore
                node_stats_line.append(len(node.proteome_ids))  # type: ignore
                node_stats.append("\t".join([str(string) for string in node_stats_line]))  # fmt:skip
                if render_tree:
                    header_f_by_node_name[node.name] = self.generate_header_for_node(node, dirs)  # fmt:skip
                charts_f_by_node_name[node.name] = self.generate_chart_for_node(node, dirs, plot_format, fontsize)  # fmt:skip
            print("[STATUS] - Writing %s ... " % node_stats_f)
            with open(node_stats_f, "w") as node_stats_fh:
                node_stats_fh.write("\n".join(node_stats) + "\n")
            print("[STATUS] - Writing %s ... " % node_clusters_f)
            with open(node_clusters_f, "w") as node_clusters_fh:
                node_clusters_fh.write("\n".join(node_clusters) + "\n")
            if render_tree:
                self.plot_tree(header_f_by_node_name, charts_f_by_node_name, dirs)
            else:
                self.plot_text_tree(dirs)
