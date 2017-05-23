#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: subset_custer_metrics_by_nodes.py        -c <FILE> -t <FILE> -a <FILE> [-g <FILE>]
                                                [-n <FLOAT>] [-d <FLOAT>] [-p <FLOAT>]
                                                [-o <STRING>]
                                                [-h|--help]

    Options:
        -h, --help                          show this
        -c, --cluster_metrics <FILE>        cluster_metrics_domains_detailed.IPR.txt file
        -t, --tree_cluster_metrics <FILE>   tree.cluster_metrics.txt file
        -a, --cluster_metrics_ALO <FILE>    Needed for getting number of proteins for clusters without domains
        -g, --node_names <FILE>             CSV file of nodes and names (only output for those nodes will be supplied)
        -n, --node_proteome_cov <FLOAT>     Minimum proteome coverage of cluster in tree.cluster_metrics.txt [default: 0.9]
                                                Clusters in output will be sorted (decreasing) by node_proteome_cov
        -d, --domain_proteome_cov <FLOAT>   Minimum proteome coverage by proteins with domain in cluster [default: 0.85]
                                                Domains of clusters in output will be sorted (decreasing) by node_proteome_cov
        -p, --domain_protein_cov <FLOAT>    Minimum protein coverage of domain in cluster [default: 0.0]
        -o, --outprefix <STRING>            Output prefix

"""

import re
import sys
import operator
from docopt import docopt
from os.path import isfile, join, exists, realpath, dirname, basename

def read_file(infile):
    if not infile or not exists(infile):
        sys.exit("[ERROR] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        for line in fh:
            line = line.replace(r'\r','\n')
            if not line.startswith("#"):
                yield line.rstrip("\n")


class DomainObj():
    def __init__(self, domain_id, domain_description, domain_proteome_cov, domain_protein_cov):
        self.domain_id = domain_id
        self.domain_description = domain_description
        self.domain_proteome_cov = domain_proteome_cov
        self.domain_protein_cov = domain_protein_cov

class ClusterObj():
    def __init__(self, cluster_id, node_id, node_description, node_label, node_proteome_cov):
        self.cluster_id = cluster_id
        self.node_proteome_cov = node_proteome_cov  # output gets sorted by node_proteome_cov
        self.node_id = node_id
        self.node_description = node_description
        self.node_label = node_label
        self.domainObjs = []
        self.protein_count = 0

    def add_domainObj(self, domainObj):
        self.domainObjs.append(domainObj)

class DataCollection():
    def __init__(self):
        self.cluster_ids_by_node_id = {}
        self.clusterObjs_by_cluster_id = {}
        self.node_description_by_node_id = {}
        self.node_label_by_node_id = {}
        self.node_ids = []

    def parse_node_names(self, node_names_f):
        print "[+] Parsing %s ..." % (node_names_f)
        for line in read_file(node_names_f):
            col = line.split(',')
            try:
                node_id = col[0]
                node_description = col[1]
                node_label = col[2]
            except IndexError:
                sys.exit("[ERROR] - File '%s' is not CSV or does not contain three columns." % (node_names_f))
            self.node_description_by_node_id[node_id] = node_description
            self.node_label_by_node_id[node_id] = node_label
            self.node_ids.append(node_id)

    def parse_tree_cluster_metrics(self, tree_cluster_metrics_f, node_proteome_cov):
        print "[+] Parsing %s ..." % (tree_cluster_metrics_f)
        for line in read_file(tree_cluster_metrics_f):
            col = line.split()
            cluster_id = col[0]
            node_id = col[1]
            try:
                node_proteome_cov = float(col[3])
                if not node_names_f:
                    if not node_id in self.node_description_by_node_id:
                        self.node_description_by_node_id[node_id] = node_id
                        self.node_ids.append(node_id)
                        self.node_label_by_node_id[node_id] = node_id
                if node_id in self.node_description_by_node_id:
                    if node_proteome_cov >= NODE_PROTEOME_COV:
                        node_description = self.node_description_by_node_id[node_id]
                        node_label = self.node_label_by_node_id[node_id]
                        clusterObj = ClusterObj(cluster_id, node_id, node_description, node_label, node_proteome_cov)
                        self.clusterObjs_by_cluster_id[cluster_id] = clusterObj
                        if node_id not in self.cluster_ids_by_node_id:
                            self.cluster_ids_by_node_id[node_id] = []
                        self.cluster_ids_by_node_id[node_id].append(cluster_id)
                        self.clusterObjs_by_cluster_id[cluster_id] = clusterObj
            except ValueError:
                pass

    def parse_cluster_metrics_ALO(self, cluster_metrics_ALO_f):
        print "[+] Parsing %s ..." % (cluster_metrics_ALO_f)
        for line in read_file(cluster_metrics_ALO_f):
            if not line.startswith("#"):
                col = line.split()
                cluster_id = col[0]
                if cluster_id in self.clusterObjs_by_cluster_id:
                    try:
                        protein_count = int(col[3])
                        self.clusterObjs_by_cluster_id[cluster_id].protein_count = protein_count
                    except ValueError:
                        pass


    def parse_cluster_metrics(self, cluster_metrics_f, DOMAIN_PROTEOME_COV, DOMAIN_PROTEIN_COV):
        print "[+] Parsing %s ..." % (cluster_metrics_f)
        for line in read_file(cluster_metrics_f):
            col = line.split("\t")
            cluster_id = col[0]
            domain_id = col[2]
            domain_description = col[3]
            protein_count = int(col[4])
            protein_count_with_domain = int(col[5])
            domain_protein_cov = protein_count/protein_count_with_domain
            domain_proteome_cov = float(col[6])
            if cluster_id in self.clusterObjs_by_cluster_id:
                if not self.clusterObjs_by_cluster_id[cluster_id].protein_count:
                    self.clusterObjs_by_cluster_id[cluster_id].protein_count = protein_count
                if domain_proteome_cov >= DOMAIN_PROTEOME_COV:
                    if domain_protein_cov >= DOMAIN_PROTEIN_COV:
                        domainObj = DomainObj(domain_id, domain_description, domain_proteome_cov, domain_protein_cov)
                        self.add_domainObj(cluster_id, domainObj)


    def add_domainObj(self, cluster_id, domainObj):
        self.clusterObjs_by_cluster_id[cluster_id].add_domainObj(domainObj)

    def write(self, prefix):
        output_cluster_domains = []
        for node_id in self.node_ids:
            if node_id in self.cluster_ids_by_node_id:
                clusterObjs = [self.clusterObjs_by_cluster_id[cluster_id] for cluster_id in self.cluster_ids_by_node_id[node_id]]
                for clusterObj in sorted(clusterObjs, key=lambda x: x.node_proteome_cov, reverse=True):
                    line_cluster_domains = []
                    line_cluster_domains.append(clusterObj.cluster_id)
                    line_cluster_domains.append(clusterObj.node_label)
                    line_cluster_domains.append(clusterObj.node_description)
                    line_cluster_domains.append(clusterObj.protein_count)
                    line_cluster_domains.append(clusterObj.node_proteome_cov)
                    if clusterObj.domainObjs:
                        domain_ids = []
                        domain_descriptions = []
                        for domainObj in sorted(clusterObj.domainObjs, key=lambda x: x.domain_proteome_cov, reverse=True):
                            domain_ids.append(domainObj.domain_id)
                            domain_descriptions.append(domainObj.domain_description)
                        line_cluster_domains.append(";".join(domain_ids))
                        line_cluster_domains.append(";".join(domain_descriptions))
                    else:
                        line_cluster_domains.append("None")
                        line_cluster_domains.append("None")
                    output_cluster_domains.append("\t".join([str(x) for x in line_cluster_domains]))
        output_f = 'cluster_domains_by_node.tsv'
        if prefix:
            output_f = prefix + "." + output_f
        print "[+] Writing %s ..." % (output_f)
        with open(output_f, 'w') as output_fh:
            output_fh.write("\n".join([str(x) for x in output_cluster_domains]))


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    node_names_f = args['--node_names']
    tree_cluster_metrics_f = args['--tree_cluster_metrics']
    cluster_metrics_f = args['--cluster_metrics']
    cluster_metrics_ALO_f = args['--cluster_metrics_ALO']
    NODE_PROTEOME_COV = float(args['--node_proteome_cov'])
    DOMAIN_PROTEOME_COV = float(args['--domain_proteome_cov'])
    DOMAIN_PROTEIN_COV = float(args['--domain_protein_cov'])
    out_prefix = args['--outprefix']

    print "[+] Start ..."
    dataCollection = DataCollection()
    if node_names_f:
        dataCollection.parse_node_names(node_names_f)
    dataCollection.parse_tree_cluster_metrics(tree_cluster_metrics_f, NODE_PROTEOME_COV)
    dataCollection.parse_cluster_metrics_ALO(cluster_metrics_ALO_f)
    dataCollection.parse_cluster_metrics(cluster_metrics_f, DOMAIN_PROTEOME_COV, DOMAIN_PROTEIN_COV)
    dataCollection.write(out_prefix)
