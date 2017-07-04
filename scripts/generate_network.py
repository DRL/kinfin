#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: generate_network.py      -m <FILE> -c <FILE> [-o <STR>] [--exclude_universal]
                                [-h|--help]

    Options:
        -h --help                               show this
        -m, --cluster_metrics <FILE>            *.cluster_metrics.txt file (e.g. TAXON.cluster_metrics.txt) from KinFin output
        -c, --config_file <FILE>                config.txt used in Kinfin analysis
        -o, --out_prefix <STR>                  Outprefix (default: graph)
        --exclude_universal                     Excludes clusters in which all taxa are present from edge weights (default: False)


"""

from __future__ import division
import sys
import networkx as nx
from os.path import basename, isfile, abspath, splitext, join, exists
from os import getcwd, mkdir
import shutil

from docopt import docopt
from itertools import combinations

def check_file(f):
    if not f or not exists(f):
        sys.exit("[-] File %s does not exist." % (f))
    else:
        return f

def read_file(f):
    with open(f) as fh:
        for line in fh:
            yield line.rstrip("\n")

def parse_cluster_stats_f(cluster_stats_f, proteomeObj_by_proteome_id, attribute, exclude_universal):
    print "[+] Parsing %s ... " % cluster_stats_f
    proteome_id_by_idx = None
    cluster_type_idx = None
    edge_weights = {}
    for line in read_file(cluster_stats_f):
        temp = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            proteome_id_by_idx = {idx: col.replace("_count", "") for idx, col in enumerate(temp) if col.replace("_count", "") in proteomeObj_by_proteome_id}
            #proteome_id_fields = [column.replace("_count", "") for column in temp if column.replace("_count", "") in proteomeObj_by_proteome_id]
            #proteome_id_idx = [idx for idx, column in enumerate(temp) if column.replace("_count", "") in proteomeObj_by_proteome_id]
            #proteome_id_by_idx = {idx : proteome_id for idx, proteome_id in zip(proteome_id_idx, sorted(proteome_id_fields))}
            if not proteome_id_by_idx:
                sys.exit("[-] No column header ending in '_count' found in %s" % (",".join(temp)))
            for idx, col in enumerate(temp):
                if col == "attribute_cluster_type":
                    cluster_type_idx = idx
            if not cluster_type_idx:
                sys.exit("[-] No column header 'cluster_type' found in %s" % (",".join(temp)))
        else:
            protein_counts_by_proteome_id = {}
            for idx, proteome_id in proteome_id_by_idx.items():
                if int(temp[idx]) > 0:
                    protein_counts_by_proteome_id[proteome_id] = int(temp[idx])
            cluster_type = temp[cluster_type_idx]
            for proteome_id, count in protein_counts_by_proteome_id.items():
                proteomeObj_by_proteome_id[proteome_id].add_cluster(cluster_type, count)
            if cluster_type  == "shared":
                # only clusters with more than one proteome
                if exclude_universal and len(protein_counts_by_proteome_id) == len(proteome_id_by_idx):
                    pass
                else:
                    for combination in combinations(protein_counts_by_proteome_id.keys(), 2):
                        edge_nodes = frozenset(combination)
                        if not edge_nodes in edge_weights:
                            edge_weights[edge_nodes] = 0
                        edge_weights[edge_nodes] += 1

    #print edge_weights
    weighted_edges = generate_edges(edge_weights)
    for proteome_id, proteomeObj in proteomeObj_by_proteome_id.items():
        proteomeObj.level_by_attribute['protein_count'] = proteomeObj.protein_counts['total']
    return weighted_edges, proteomeObj_by_proteome_id

def generate_edges(edge_weights):
    weighted_edges = []
    max_edge_weight = max(edge_weights.values())
    print "[+] Max edge weight is %s, ..." % (max_edge_weight)
    for edge_nodes, weight in edge_weights.items():
        nodes = list(edge_nodes)
        edge_tuple = (nodes[0], nodes[1], weight)
        weighted_edges.append(edge_tuple)
    return weighted_edges

def generate_outpath_by_name(out_prefix):
    outpath_by_name = {}
    string = ''
    if out_prefix:
        string = "%s.graph" % out_prefix
    else:
        string = "graph"
    outpath_by_name["graphml_f"] = "%s.graphml" % (string)
    outpath_by_name["gex_f"] =  "%s.gexf" % (string)
    return outpath_by_name


def construct_graphs(weighted_edges, proteomeObj_by_proteome_id, attributes):
    print "[+] Building graphs"
    G = nx.Graph()
    G.name = "Graph"
    proteome_id_by_idx = {}
    for proteome_id, proteomeObj  in proteomeObj_by_proteome_id.items():
        G.add_node(proteomeObj.proteome_id, {k :v for k, v in proteomeObj.level_by_attribute.items()})
        G.add_node(proteomeObj.proteome_id, {k :v for k, v in proteomeObj.level_by_attribute.items()})
    G.add_weighted_edges_from(weighted_edges)
    print nx.info(G)
    print "[+] Saving network %s" % outpath_by_name["graphml_f"]
    nx.write_graphml(G, outpath_by_name["graphml_f"])
    print "[+] Saving network %s" % outpath_by_name["gex_f"]
    nx.write_gexf(G, outpath_by_name["gex_f"])


def santisise_args(args):
    sane_args = {}
    sane_args['--cluster_metrics'] = check_file(args['--cluster_metrics'])
    sane_args['--config_file'] = check_file(args['--config_file'])
    sane_args['--out_prefix'] = args['--out_prefix']
    return sane_args

def parse_classification_f(species_classification_f):
    print "[+] Parsing SpeciesClassification file: %s ..." % (species_classification_f)
    attributes = []
    proteomeObj_by_proteome_id = {}
    for line in read_file(species_classification_f):
        if line.startswith("#"):
            if not attributes:
                attributes = [x.strip() for x in line.lstrip("#").split(",")]
                if not 'IDX' == attributes[0] or not 'TAXON' == attributes[1]:
                    sys.exit("[-] First/second element have to be IDX/TAXON.\n\t%s" % (attributes))
                else:
                    pass # accounts for SpeciesIDs that are commented out for Orthifinder
        elif line.strip():
            temp = line.split(",")
            if not len(temp) == len(attributes):
                sys.exit("[-] number of columns in line differs from header\n\t%s\n\t%s" % (attributes, temp))
            if temp[1] in proteomeObj_by_proteome_id:
                sys.exit("[-] 'proteome' should be unique. %s was encountered multiple times" % (temp[0]))
            proteome_idx = temp[0]
            proteome_id = temp[1]
            level_by_attribute = {x : '' for x in attributes}
            for idx, level in enumerate(temp):
                attribute = attributes[idx]
                level_by_attribute[attribute] = level
            proteomeObj = ProteomeObj(proteome_id, proteome_idx, level_by_attribute)
            proteomeObj_by_proteome_id[proteome_id] = proteomeObj
        else:
            pass
    return proteomeObj_by_proteome_id, attributes


class ProteomeObj():
    def __init__(self, proteome_id, proteome_idx, level_by_attribute):
        self.proteome_idx = proteome_idx
        self.proteome_id = proteome_id
        self.level_by_attribute = level_by_attribute

        self.protein_counts = {"total": 0, "singleton" : 0, "shared" : 0, "specific" : 0}
        self.cluster_counts = {"total": 0, "singleton" : 0, "shared" : 0, "specific" : 0}

    def add_cluster(self, cluster_type, count):
        self.protein_counts[cluster_type] += count
        self.protein_counts["total"] += count
        self.cluster_counts[cluster_type] += 1
        self.cluster_counts["total"] += 1

if __name__ == "__main__":
    __version__ = 0.3

    args = docopt(__doc__)
    exclude_universal = args['--exclude_universal']
    sane_args = santisise_args(args)
    cluster_stats_f = sane_args['--cluster_metrics']
    species_classification_f = sane_args['--config_file']
    out_prefix = sane_args['--out_prefix']

    #######################################################################

    proteomeObj_by_proteome_id, attributes = parse_classification_f(species_classification_f)
    outpath_by_name = generate_outpath_by_name(out_prefix)
    weighted_edges, proteomeObj = parse_cluster_stats_f(cluster_stats_f, proteomeObj_by_proteome_id, attributes, exclude_universal)
    construct_graphs(weighted_edges, proteomeObj_by_proteome_id, attributes)
