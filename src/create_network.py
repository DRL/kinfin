#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: generate_and_analyse_network.py      -c <FILE> -s <FILE> -a <STR>
                                            [-n <STR>] [-t <FLOAT>] [-e <STR>]
                                            [--colour <FILE>]
                                            [-o <STR>] [-d]
                                            [-h|--help]

    Options:
        -h --help                               show this
        -c, --cluster_stats <FILE>              Kinfin cluster_stats.txt output (i.e. PROTEOME.cluster_stats.txt)
        -s, --species_classification <FILE>     SpeciesClassification.txt used in Kinfin analysis
        -a, --attribute <STR>                   Attribute to use as label (colouring and legend)
        -n, --normalisation <STR>               Procedure to normalise weights of edges (default: none)
                                                    - "max_weight"
        -t, --weight_threshold <FLOAT>          Minimal weight threshold (default: 0.0)
        -e, --exclude_proteome_ids <STR>        Proteomes to be excluded (seprated by ",")
        -d, --draw_all_edges                    Draw all edges (as opposed to only filtered ones)
        --colour_f <FILE>                       CSV file containing colours of attribute levels
        -o, --out_prefix <STR>                  Outprefix (default: graph)

"""

from __future__ import division
import sys
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from networkx.algorithms import tree
from os.path import basename, isfile, abspath, splitext, join, exists
from os import getcwd, mkdir
import shutil

from docopt import docopt
from itertools import combinations
from math import log
import scipy
import numpy as np
import matplotlib as mat
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter, NullFormatter
from mpl_toolkits.mplot3d import Axes3D

mat.use('agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import seaborn as sns
sns.set_context("talk")
sns.set(style="whitegrid")
import pylab
sns.set_color_codes("pastel")
mat.rc('ytick', labelsize=20)
mat.rc('xtick', labelsize=20)
axis_font = {'size':'20'}
mat.rcParams.update({'font.size': 22})

def check_file(f):
    if not f or not exists(f):
        sys.exit("[ERROR] - File %s does not exist." % (f))
    else:
        return f

def read_file(f):
    with open(f) as fh:
        for line in fh:
            yield line.rstrip("\n")

class ProteomeObj():
    def __init__(self, proteome_id, label, colour, idx):
        self.proteome_id = proteome_id
        self.label = label
        self.colour = colour
        self.protein_counts = {"total": 0, "singleton" : 0, "shared" : 0, "specific" : 0}
        self.cluster_counts = {"total": 0, "singleton" : 0, "shared" : 0, "specific" : 0}
        self.idx = idx

    def add_cluster(self, cluster_type, count):
        self.protein_counts[cluster_type] += count
        self.protein_counts["total"] += count
        self.cluster_counts[cluster_type] += 1
        self.cluster_counts["total"] += 1


def parse_species_classification(species_classification_f):
    label_by_proteome_id = {}
    proteome_ids_by_label = {}
    target_attribute_idx = None
    proteome_ids = []
    for line in read_file(species_classification_f):
        if line.startswith("#"):
            attributes = [x.strip() for x in line.lstrip("#").split(",")]
            if not target_attribute in set(attributes):
                sys.exit("[ERROR] - Attribute %s not found in %s" % (target_attribute, ",".join(attributes)))
            target_attribute_idx = attributes.index(target_attribute)
        elif line.strip():
            temp = line.split(",")
            if not len(temp) == len(attributes):
                sys.exit("[ERROR] - number of columns in line differs from header\n\t%s\n\t%s" % (attributes, temp))
            proteome_id = temp[1]
            proteome_ids.append(proteome_id)
            label = temp[target_attribute_idx]
            label_by_proteome_id[proteome_id] = label
            if not label in proteome_ids_by_label:
                proteome_ids_by_label[label] = []
            proteome_ids_by_label[label].append(proteome_id)
        else:
            pass
    return proteome_ids, label_by_proteome_id, proteome_ids_by_label

def parse_cluster_stats_f(cluster_stats_f):
    print "[STATUS] - Parsing %s ... " % cluster_stats_f
    proteome_id_by_idx = None
    cluster_type_idx = None
    edge_count = {}
    cluster_count = {}
    weight_by_all_edges = {} # initiated once proteomeIDs are known
    weight_by_positional_edges = {} # initiated on demand
    for line in read_file(cluster_stats_f):
        temp = line.rstrip("\n").split()
        if line.startswith("#"):
            proteome_id_by_idx = {idx: col.replace("_count", "") for idx, col in enumerate(temp) if col.replace("_count", "") in proteome_ids}
            if not proteome_id_by_idx:
                sys.exit("[ERROR] - No column header ending in '_count' found in %s" % (",".join(temp)))
            for idx, col in enumerate(temp):
                if col == "attribute_cluster_type":
                    cluster_type_idx = idx
            if not cluster_type_idx:
                sys.exit("[ERROR] - No column header 'cluster_type' found in %s" % (",".join(temp)))

            # initiate weight_by_all_edges
            for combination in combinations(proteome_id_by_idx.values(), 2):
                edge = frozenset(combination)
                weight_by_all_edges[edge] = 0
        else:
            protein_counts_by_proteome_id = {}
            for idx, proteome_id in proteome_id_by_idx.items():
                if int(temp[idx]) > 0:
                    protein_counts_by_proteome_id[proteome_id] = int(temp[idx])
            cluster_type = temp[cluster_type_idx]
            proteome_ids_present = set(protein_counts_by_proteome_id.keys())
            for proteome_id, count in protein_counts_by_proteome_id.items():
                proteomeObj_by_proteome_id[proteome_id].add_cluster(cluster_type, count)
            if cluster_type == "shared":
                cluster_count['total'] = cluster_count.get('total', 0) + 1
                overlap_with_excluded = exclude_proteome_ids.intersection(proteome_ids_present)
                # only clusters with more than one proteome
                for combination in combinations(proteome_ids_present, 2):
                    edge = frozenset(combination)
                    edge_count['total'] = edge_count.get('total', 0) + 1
                    if not overlap_with_excluded:
                        weight_by_positional_edges[edge] = weight_by_positional_edges.get(edge, 0) + 1
                    else:
                        edge_count['excluded_proteome'] = edge_count.get('excluded_proteome', 0) + 1
                    weight_by_all_edges[edge] = weight_by_all_edges.get(edge, 0) + 1
                if overlap_with_excluded:
                    cluster_count['excluded_proteome'] = cluster_count.get('excluded_proteome', 0) + 1
    edges_all = generate_edges("none", weight_by_all_edges)
    edges_positional = generate_edges(normalisation, weight_by_positional_edges)
    return edges_all, edges_positional, cluster_count, edge_count

def generate_edges(normalisation, weight_by_positional_edges):
    edges_positional = []
    if normalisation == "none":
        for edge, weight in weight_by_positional_edges.items():
            nodes = list(edge)
            weight = weight_by_positional_edges[edge]
            edge_tuple = (nodes[0],nodes[1], weight)
            edges_positional.append(edge_tuple)
    elif normalisation == "max_weight":
        max_weight_positional = max(weight_by_positional_edges.values())
        for edge, weight in weight_by_positional_edges.items():
            nodes = list(edge)
            normalised_weight = weight/max_weight_positional
            edge_tuple = (nodes[0],nodes[1], normalised_weight)
            edges_positional.append(edge_tuple)
    else:
        pass
    return edges_positional


def generate_proteomeObjs():
    proteomeObj_by_proteome_id = {}
    labels = [label for label in proteome_ids_by_label.keys()]
    colour_by_label = {label : colors.rgb2hex(plt.cm.rainbow(idx/len(labels))) for idx, label in enumerate(labels)}
    for idx, proteome_id in enumerate(proteome_ids):
        label = label_by_proteome_id[proteome_id]
        colour = colour_by_label[label]
        proteomeObj = ProteomeObj(proteome_id, label, colour, idx)
        proteomeObj_by_proteome_id[proteome_id] = proteomeObj
    return proteomeObj_by_proteome_id

def check_normalisation(normalisation_algo):
    if not normalisation_algo in NORMALISATION_ALGOS:
        sys.exit("[ERROR] - Normalisation algorithm %s not found in %s" % (normalisation_algo, ",".join(list(NORMALISATION_ALGOS))))
    else:
        return normalisation_algo

def check_exclude_proteome_ids(exclude_proteome_ids):
    return set([proteome_id for proteome_id in exclude_proteome_ids.split(",")])


def generate_outpath_by_name(target_attribute, normalisation, weight_threshold, out_prefix):
    outpath_by_name = {}
    string = ''
    if out_prefix:
        string = "%s." % out_prefix
    outpath_by_name["main_dir"] = join(getcwd(), "%skinfin_network_analysis" % (string))
    print "[STATUS] - Output directory is \n\t%s" % (outpath_by_name["main_dir"])
    #if exists(outpath_by_name["main_dir"]):
    #    print "[STATUS] - Directory exists. Deleting directory ..."
    #    shutil.rmtree(outpath_by_name["main_dir"])
    #print "[STATUS] - Creating directories ..."
    #mkdir(outpath_by_name["main_dir"])
    outpath_by_name["nodes_dir"] = join(outpath_by_name["main_dir"], "nodes")
    if not exists(outpath_by_name["nodes_dir"]):
        print "\t%s" % (outpath_by_name["nodes_dir"])
        mkdir(outpath_by_name["nodes_dir"])
    if normalisation == "none":
        outpath_by_name["graph_f"] = join(outpath_by_name["main_dir"], "graph.pdf")
        outpath_by_name["gml_f"] = join(outpath_by_name["main_dir"], "graph.gml")
        outpath_by_name["edge_list_f"] = join(outpath_by_name["main_dir"], "graph.edge_list.txt")
    else:
        outpath_by_name["graph_f"] = join(outpath_by_name["main_dir"], "graph.%s.%s.pdf" % (normalisation, weight_threshold))
        outpath_by_name["gml_f"] = join(outpath_by_name["main_dir"], "graph.%s.%s.gml" % (normalisation, weight_threshold))
        outpath_by_name["edge_list_f"] = join(outpath_by_name["main_dir"], "graphedge_list.%s.%s.txt" % (normalisation, weight_threshold))
    outpath_by_name["plot_node_degree_f"] = join(outpath_by_name["main_dir"], "node_degree.pdf")
    return outpath_by_name

def mean(lst):
    if lst:
        return float(sum(lst)) / len(lst)
    else:
        return 0.0


def construct_graphs(edges_all, edges_positional, edges_positional_filtered, draw_all_edges):
    print "[STATUS] - Building graphs"
    G_edges_all = nx.Graph()
    G_edges_all.name = "edges_all"
    G_edges_positional = nx.Graph()
    G_edges_positional.name = "edges_positional"
    G_edges_positional_filtered = nx.Graph()
    G_edges_positional_filtered.name = "edges_positional_filtered"
    proteome_id_by_idx = {}
    for idx, proteome_id in enumerate(proteome_ids):
        #print idx, proteome_id
        # populate proteome_id_by_idx
        proteomeObj = proteomeObj_by_proteome_id[proteome_id]
        proteome_id_by_idx[proteome_id] = idx
        # add nodes
        G_edges_all.add_node(proteomeObj.proteome_id, {k :v for k, v in proteomeObj.__dict__.items()})
        G_edges_positional.add_node(proteomeObj.proteome_id, {k :v for k, v in proteomeObj.__dict__.items()})
        G_edges_positional_filtered.add_node(proteomeObj.proteome_id, {k :v for k, v in proteomeObj.__dict__.items()})
    # set max_node size
    max_node_size = max([G_edges_all.node[node]['protein_counts']['total'] for node in G_edges_all.nodes()])
    # set node sizes scaled by max_node_size*1000
    node_sizes_of_nodes = {node : 1000*G_edges_all.node[node]['protein_counts']['total']/max_node_size for node in G_edges_all.nodes()}
    # set node colour
    node_colour_of_nodes = {node : G_edges_all.node[node]['colour'] for node in G_edges_all.nodes()}
    # set node label
    node_label_of_nodes = {node : G_edges_all.node[node]['label'] for node in G_edges_all.nodes()}

    # Add edges to G_edges_all
    G_edges_all.add_weighted_edges_from(edges_all)
    # calculate node_weighted_degrees_of_nodes from G_edges_all
    node_weighted_degrees_of_nodes_all = nx.degree(G_edges_all, weight='weight')
    print "edges in edges_all w/o connections"
    for edge in edges_all:
        if edge[2] == 0:
            print edge
    print nx.info(G_edges_all)
    # Add edges to G_edges_positional
    G_edges_positional.add_weighted_edges_from(edges_positional)
    # calculate node_weighted_degrees_of_nodes from G_edges_all
    node_weighted_degrees_of_nodes_positional = nx.degree(G_edges_positional, weight='weight')
    print nx.info(G_edges_positional)

    # Add edges to G_edges_positional_filtered
    G_edges_positional_filtered.add_weighted_edges_from(edges_positional_filtered)
    # calculate node_weighted_degrees_of_nodes from G_edges_positional_filtered
    node_weighted_degrees_of_nodes_positional = nx.degree(G_edges_positional_filtered, weight='weight')
    print nx.info(G_edges_positional_filtered)

    # infer positions of nodes
    node_position_of_nodes_positional_filtered = graphviz_layout(G_edges_positional_filtered, prog="neato")
    # draw network
    draw_graph(G_edges_all, G_edges_positional, G_edges_positional_filtered, draw_all_edges, node_position_of_nodes_positional_filtered, node_sizes_of_nodes, node_colour_of_nodes, node_label_of_nodes, proteome_id_by_idx)

def draw_graph(G_edges_all, G_edges_positional, G_edges_positional_filtered, draw_all_edges, node_position_of_nodes_positional_filtered, node_sizes_of_nodes, node_colour_of_nodes, node_label_of_nodes, proteome_id_by_idx):
    f, ax = plt.subplots(figsize=FIGSIZE)
    for label, proteome_ids in proteome_ids_by_label.items():
        node_list = proteome_ids
        node_sizes = [node_sizes_of_nodes[node] for node in node_list]
        node_colour = [node_colour_of_nodes[node] for node in node_list]
        node_label = label
        nx.draw_networkx_nodes(G_edges_all, node_position_of_nodes_positional_filtered, nodelist=node_list, node_size=node_sizes, node_color=node_colour, label=node_label)

    G_edges_blessed = []
    G_edges_to_excluded = []
    G_edges_below_threshold = []
    for edge in G_edges_all.edges(data='weight'):
        if not nodes_connected(G_edges_positional, edge[0], edge[1]):
            G_edges_to_excluded.append(edge)
        else:
            if not nodes_connected(G_edges_positional_filtered, edge[0], edge[1]):
                G_edges_below_threshold.append(edge)
            else:
                G_edges_blessed.append(edge)
    #edges_all, weights_all = zip(*nx.get_edge_attributes(G_edges_all, 'weight').items())
    #edges_positional_filtered, weights_positional_filtered = zip(*nx.get_edge_attributes(G_edges_positional, 'weight').items())
    #edges_all_labels = nx.get_edge_attributes(G_edges_all, 'weight')
    #edges_positional_labels = nx.get_edge_attributes(G_edges_positional, 'weight')
    #edges_positional_filtered_labels = nx.get_edge_attributes(G_edges_positional_filtered, 'weight')
    #nx.draw_networkx_edge_labels(G_edges_all, node_position_of_nodes_positional_filtered, edge_labels=edges_positional_labels)

    ## With weighted edges
    #edges_positional, weights_positional = zip(*nx.get_edge_attributes(G_edges_positional, 'weight').items())
    #edges, weights = zip(*nx.get_edge_attributes(G_edges_all, 'weight').items())
    #nx_edges = nx.draw_networkx_edges(G_edges_positional_filtered, node_position_of_nodes_positional_filtered, edgelist=edges, norm=colors.LogNorm(vmin=np.array(weights).min(), vmax=np.array(weights).max()), edge_color=weights, alpha=0.8, width=2, edge_cmap=plt.cm.Greys)
    #nx.draw_networkx_labels(G_edges_all, node_position_of_nodes_positional_filtered, labels=proteome_id_by_idx)

    ## With three edge colours
    nx.draw_networkx_edges(G_edges_all, node_position_of_nodes_positional_filtered, edgelist=G_edges_to_excluded, edge_color="lightblue", alpha=0.3, width=1)
    nx.draw_networkx_edges(G_edges_positional, node_position_of_nodes_positional_filtered, edgelist=G_edges_below_threshold, edge_color="lightgrey", alpha=0.3, width=1)
    nx.draw_networkx_edges(G_edges_positional_filtered, node_position_of_nodes_positional_filtered, edgelist=G_edges_blessed, edge_color="black", alpha=0.3, width=1)
    nx.draw_networkx_labels(G_edges_all, node_position_of_nodes_positional_filtered, labels=proteome_id_by_idx)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), labelspacing=1, borderaxespad=0., mode="expand", scatterpoints=1)
    for legendHandle in legend.legendHandles:
        legendHandle._sizes = [100]

    #plt.colorbar(nx_edges, orientation="horizontal")
    plt.axis('off')
    f.savefig(outpath_by_name["graph_f"], dpi=1000)
    print "[STATUS] - Saving network %s" % outpath_by_name["graph_f"]

    plt.close()
    sys.exit()


    f, ax = plt.subplots(figsize=FIGSIZE)
    for proteome_id, p_rank in sorted(pagerank.items(), key=lambda x: x[1], reverse=True):
        print proteome_id, p_rank
    node_size = {}
    max_proteome_size = max([proteomeObj_by_proteome_id[proteome_id].protein_counts['total'] for proteome_id in proteomeObj_by_proteome_id])
    for proteomeObj in proteomeObjs:
        superphylum = proteomeObj.superphylum
        if not superphylum in superphyla_d:
            superphyla_d[superphylum] = []
        if not superphylum in node_size:
            node_size[superphylum] = []
        superphyla_d[superphylum].append(proteomeObj.proteome_id)
        node_size[superphylum].append(1000*proteomeObj.protein_counts['total']/max_proteome_size)
        colour_d[superphylum] = colour_by_superphylum[superphylum]


    for node, degree in sorted(all_edges_G.degree(weight='weight').items(), key=lambda x: x[1], reverse=True):
        print node, degree
    for superphylum in superphyla_d:
        #print superphylum, superphyla_d[superphylum]
        nx.draw_networkx_nodes(some_edges_G, pos, nodelist=superphyla_d[superphylum], node_size=node_size[superphylum], node_color=colour_d[superphylum], label=superphylum)

    print "edges to be drawn", len(edgelist)
    #print edgelist
    nx.draw_networkx_edges(some_edges_G, pos, edge_color="darkgrey", alpha=0.1, width=1)
    print nx.info(all_edges_G)
    print nx.info(some_edges_G)

    nx.draw_networkx_labels(some_edges_G, pos, labels)
    #plt.colorbar(edges)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), labelspacing=1, borderaxespad=0., mode="expand", scatterpoints=1)
    for legendHandle in legend.legendHandles:
        legendHandle._sizes = [100]
    # for MCL
    #nx.write_edgelist(all_edges_G,'50helminths.adjacency_list.txt', data=['weight'])
    plt.axis('off')
    nx_f = graph_f
    f.savefig(nx_f, dpi=1000)
    print "[STATUS] - Saving network %s" % nx_f
    nx.write_gml(all_edges_G,'50helminths.adjacency_list.gml')
    nx.write_graphml(all_edges_G,'50helminths.adjacency_list.xml')
    plt.close()

def nodes_connected(G, u, v):
    return u in G.neighbors(v)

def filter_edges(edges_positional, weight_threshold, edge_count):
    edges_positional_filtered = []
    for edge in edges_positional:
        if edge[2] >= weight_threshold:
            edges_positional_filtered.append(edge)
        else:
             edge_count['excluded_weight'] = edge_count.get('excluded_weight', 0) + 1
    return edges_positional_filtered, edge_count

def santisise_args(args):
    sane_args = {}
    sane_args['--cluster_stats'] = check_file(args['--cluster_stats'])
    sane_args['--colour_f'] = check_file(args['--colour_f'])
    sane_args['--species_classification'] = check_file(args['--species_classification'])
    sane_args['--attribute'] = args['--attribute']
    sane_args['--normalisation'] = check_normalisation(args['--normalisation'])
    sane_args['--weight_threshold'] = float(args['--weight_threshold'])
    if args['--exclude_proteome_ids']:
        sane_args['--exclude_proteome_ids'] = check_exclude_proteome_ids(args['--exclude_proteome_ids'])
    else:
        sane_args['--exclude_proteome_ids'] = set()
    sane_args['--draw_all_edges'] = args['--draw_all_edges']
    sane_args['--out_prefix'] = args['--out_prefix']
    return sane_args

if __name__ == "__main__":
    '''
    IDEAS:
        - get pos for nodes with filtered edges, but draw with unfiltered edges
        - draw nodes as pie charts
        - Report on number of connections for each proteome
            - draw degree distribution
        - Filter proteomes to exclude in filter edges
        - that way they get drawn when -d
         src/create_network.py -c ~/Dropbox/projects/manuscripts/50helminths/results/50helminth.20161115.master.kinfin_results/PROTEOME/PROTEOME.cluster_stats.txt -s ~/Dropbox/projects/manuscripts/50helminths/results/50helminth.20161115.master.kinfin_results/PROTEOME/PROTEOME.cluster_stats.txt -a superphylum -n "max_weight" -t 0.30 -o test -d

    '''
    __version__ = 0.2

    NORMALISATION_ALGOS = set(["max_weight", "none"])
    FIGSIZE = (16,16)

    args = docopt(__doc__)
    sane_args = santisise_args(args)
    cluster_stats_f = sane_args['--cluster_stats']
    species_classification_f = sane_args['--species_classification']
    target_attribute = sane_args['--attribute']
    normalisation = sane_args['--normalisation']
    weight_threshold = sane_args['--weight_threshold']
    exclude_proteome_ids = sane_args['--exclude_proteome_ids']
    draw_all_edges = sane_args['--draw_all_edges']
    out_prefix = sane_args['--out_prefix']
    colour_f = sane_args['--colour_f']

    results = {}
    outpath_by_name = generate_outpath_by_name(target_attribute, normalisation, weight_threshold, out_prefix)

    proteome_ids, label_by_proteome_id, proteome_ids_by_label = parse_species_classification(species_classification_f)
    colour_by_proteome_id = parse_colour_f(colour_f)
    results['proteome_count'] = len(proteome_ids)
    results['labels'] = ",".join(["%s (%s)" % (label, len(proteomes)) for label, proteomes in proteome_ids_by_label.items()])
    proteomeObj_by_proteome_id = generate_proteomeObjs()

    #print proteome_ids
    edges_all, edges_positional, cluster_count, edge_count = parse_cluster_stats_f(cluster_stats_f)
    print "[RESULT] - Edges excluded"
    edges_positional_filtered, edge_count = filter_edges(edges_positional, weight_threshold, edge_count)

    construct_graphs(edges_all, edges_positional, edges_positional_filtered, draw_all_edges)
