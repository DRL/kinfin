#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: filter_goatools_enrichment.py        -i <FILE> -t <FILE> -n <STR> -g <FILE> [-o <STR>]
                                            [-m <FLOAT>] [--BP] [--CC]
                                            [--print_depleted]
                                            [-h|--help]

    Options:
        -h --help                               show this
        -i, --infile <FILE>                     Goatools enrichment file (determines proteome membership by prefix)
        -g, --groups <FILE>                     Orthologous groups files
        -t, --node_metrics <FILE>               KinFin's tree.node_metrics.txt
        -n, --node_id <STR>                     Node ID
        -o, --out_prefix <STR>                  Outprefix
        -m, --min_proteome_coverage <FLOAT>     Minimum percentage of proteomes with GO term [default: 0.0]
        --BP                                    Print BPs
        --CC                                    Print CCs
        --print_depleted                        Print depleted [default: False]
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

def parse_node_metrics(node_metrics_f):
    proteome_count_by_node_id = {}
    print "[+] Parsing Node-Metrics file: %s ..." % (node_metrics_f)
    for line in read_file(node_metrics_f):
        if line.startswith("#") or line.startswith("nodeID"):
            pass
        elif line.strip():
            temp = line.split()
            node_id = temp[0]
            proteome_count = int(temp[-1])
            proteome_count_by_node_id[node_id] = proteome_count
        else:
            pass
    return proteome_count_by_node_id

def parse_groups(groups_f):
    cluster_id_by_protein_id = {}
    cluster_size_by_cluster_id = {}
    print "[+] Parsing Groups file: %s ..." % (groups_f)
    for line in read_file(groups_f):
        temp = line.split()
        cluster_id = temp[0].replace(":", "")
        protein_ids = temp[1:]
        cluster_size_by_cluster_id[cluster_id] = len(protein_ids)
        for protein_id in protein_ids:
            cluster_id_by_protein_id[protein_id] = cluster_id
    return cluster_id_by_protein_id, cluster_size_by_cluster_id

def parse_enrichment(enrichment_f):
    output = []
    print "[+] Parsing Enrichment file: %s ..." % (enrichment_f)
    for line in read_file(enrichment_f):
        if line.startswith("#"):
            temp = line.split("\t")
            temp[-1] = "proteome_coverage"
            temp.append("proteomes")
            temp.append("clusters")
            output_line = "\t".join(temp)
            output.append(output_line)
        else:
            temp = line.split("\t")
            NS = temp[1]
            enrichment = temp[2]
            if NS in NS_parsable and enrichment in enrichment_parsable:
                if temp[10]:
                    protein_ids = [protein_id for protein_id in temp[10].split(", ")]
                    proteome_ids = set([protein_id.split("|")[0] for protein_id in protein_ids])
                    proteome_coverage = len(proteome_ids)/proteome_count_by_node_id[node_id]
                    if proteome_coverage >= min_coverage:
                        temp[10] = str(proteome_coverage)
                        temp.append("%s/%s" % (len(proteome_ids), proteome_count_by_node_id[node_id]))
                        protein_ids_by_cluster_id = {}
                        for protein_id in protein_ids:
                            cluster_id = cluster_id_by_protein_id[protein_id]
                            if not cluster_id in protein_ids_by_cluster_id:
                                protein_ids_by_cluster_id[cluster_id] = []
                            protein_ids_by_cluster_id[cluster_id].append(protein_id)
                        cluster_string = ";".join(["%s:%s" % (cluster_id, len(protein_ids_by_cluster_id[cluster_id])/cluster_size_by_cluster_id[cluster_id]) for cluster_id in sorted(protein_ids_by_cluster_id)])
                        temp.append(cluster_string)
                else:
                    temp[10] = str(0.0)
                    temp.append("0/%s" % (proteome_count_by_node_id[node_id]))
                    temp.append("N/A")
                output_line = "\t".join(temp)
                output.append(output_line)
    output_temp = enrichment_f.split(".")[0:-1]
    output_temp.append("filtered")
    output_f = ''
    if out_prefix:
        output_f = out_prefix
    output_f += ".".join(output_temp) + ".tsv"
    with open(output_f, 'w') as output_fh:
        output_fh.write("\n".join(output) + "\n")

def santisise_args(args):
    sane_args = {}
    sane_args['--infile'] = check_file(args['--infile'])
    sane_args['--node_metrics'] = check_file(args['--node_metrics'])
    sane_args['--out_prefix'] = args['--out_prefix']
    sane_args['--min_proteome_coverage'] = float(args['--min_proteome_coverage'])
    sane_args['--node_id'] = args['--node_id']
    sane_args['--groups'] = args['--groups']
    sane_args['NS_parsable'] = set(['MF'])
    if args['--BP']:
        sane_args['NS_parsable'].add("BP")
    if args['--CC']:
        sane_args['NS_parsable'].add("CC")
    sane_args['enrichment'] = set(['e'])
    if args['--print_depleted']:
        sane_args['enrichment'].add('p')
    return sane_args

if __name__ == "__main__":
    __version__ = 0.1

    args = docopt(__doc__)
    sane_args = santisise_args(args)
    enrichment_f = sane_args['--infile']
    node_metrics_f = sane_args['--node_metrics']
    groups_f = sane_args['--groups']
    node_id = sane_args['--node_id']
    out_prefix = sane_args['--out_prefix']
    min_coverage = sane_args['--min_proteome_coverage']
    NS_parsable = sane_args['NS_parsable']
    enrichment_parsable = sane_args['enrichment']

    #######################################################################

    cluster_id_by_protein_id, cluster_size_by_cluster_id = parse_groups(groups_f)
    proteome_count_by_node_id = parse_node_metrics(node_metrics_f)
    parse_enrichment(enrichment_f)

