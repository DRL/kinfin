#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: get_count_matrix.py           -c <FILE> -t <FILE> -m <FILE> -i <FILE>
                                        [-o <STR>] [-p <FLOAT>]
                                        [-h|--help]

    Options:
        -h --help                           show this
        -o, --outprefix <STRING>            Output prefix [default: results]
        -i, --cluster_ids <FILE>            List of cluster ids
        -t, --taxon_order <FILE>            Taxon order file
        -m, --cluster_metrics <FILE>        Cluster metrics file
        -c, --clustering <FILE>             clustering file
        -p, --domain_proteome_cov <FLOAT>   Proteome coverage of domain [default: 0.75]
"""

from __future__ import division
import re
import sys
import operator
from docopt import docopt
from collections import Counter
import os

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    print "[+] Parsing file %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            if not line.startswith("#"):
                yield line.rstrip("\n")


class DomainObj():
    def __init__(self, domain_id, domain_description, domain_proteome_cov, domain_protein_cov):
        self.domain_id = domain_id
        self.domain_description = domain_description
        self.domain_proteome_cov = domain_proteome_cov
        self.domain_protein_cov = domain_protein_cov


class ClusterObj():
    def __init__(self, cluster_id, protein_ids):
        self.cluster_id = cluster_id
        self.protein_ids = set(protein_ids)
        self.protein_count = len(protein_ids)
        self.protein_count_by_taxon_id = self.protein_count_by_taxon_id(protein_ids)
        self.taxon_count = len(self.protein_count_by_taxon_id)
        self.cluster_type = self.get_cluster_type()
        self.domainObjs = []

    def protein_count_by_taxon_id(self, protein_ids):
        taxon_ids = [protein_id.split(".")[0] for protein_id in protein_ids]
        return Counter(taxon_ids)

    def add_domainObj(self, domainObj):
        self.domainObjs.append(domainObj)

    def get_cluster_type(self):
        if self.taxon_count > 1:
            return "multitaxon"
        elif self.protein_count > 1:
            return "monotaxon"
        elif self.protein_count == 1:
            return "singleton"
        else:
            return None


class DataCollection():
    def __init__(self, outprefix):
        self.outprefix = outprefix
        self.taxon_ids = []
        self.cluster_ids = []
        self.clusterObj_by_cluster_id = {}


    def parse_taxon_order_f(self, taxon_order_f):
        for line in read_file(taxon_order_f):
            self.taxon_ids.append(line)

    def parse_target_cluster_id_f(self, target_cluster_id_f):
        for line in read_file(target_cluster_id_f):
            cluster_id = line
            self.cluster_ids.append(cluster_id)
            self.clusterObj_by_cluster_id[cluster_id] = {}

    def parse_clustering_f(self, clustering_f):
        for line in read_file(clustering_f):
            cluster_id, protein_string = line.rstrip("\n").split(": ")
            protein_ids = protein_string.split(" ")
            if cluster_id in self.clusterObj_by_cluster_id:
                clusterObj = ClusterObj(cluster_id, protein_ids)
                self.clusterObj_by_cluster_id[cluster_id] = clusterObj

    def parse_cluster_metrics(self, cluster_metrics_f, DOMAIN_PROTEOME_COV):
        for line in read_file(cluster_metrics_f):
            col = line.split("\t")
            cluster_id = col[0]
            domain_id = col[2]
            domain_description = col[3].replace(",", "").replace(";", "")
            protein_count = int(col[4])
            protein_count_with_domain = int(col[5])
            domain_protein_cov = protein_count_with_domain / protein_count
            domain_proteome_cov = float(col[6])
            if cluster_id in self.clusterObj_by_cluster_id:
                if domain_proteome_cov >= DOMAIN_PROTEOME_COV:
                    domainObj = DomainObj(domain_id, domain_description, domain_proteome_cov, domain_protein_cov)
                    self.add_domainObj(cluster_id, domainObj)

    def add_domainObj(self, cluster_id, domainObj):
        self.clusterObj_by_cluster_id[cluster_id].add_domainObj(domainObj)

    def write_output(self):
        output = []
        matrix_f = "matrix.txt"
        if outprefix:
            matrix_f = "%s.%s" % (outprefix, matrix_f)
        output.append("cluster_id,protein_count,%s,Domains" % ",".join(self.taxon_ids))
        for cluster_id in self.cluster_ids:
            outline = []
            cluster_id = cluster_id
            clusterObj = self.clusterObj_by_cluster_id[cluster_id]
            outline.append(cluster_id)
            outline.append(clusterObj.protein_count)
            for taxon_id in self.taxon_ids:
                taxon_count = clusterObj.protein_count_by_taxon_id.get(taxon_id, 0)
                outline.append(taxon_count)
            domain_string = []
            if cluster_id == "OG0001332":
                for domainObj in clusterObj.domainObjs:
                    print domainObj.__dict__
            if clusterObj.domainObjs:
                for domainObj in sorted(clusterObj.domainObjs, key=lambda x: x.domain_protein_cov, reverse=True):
                    domain_string.append("%s=%s (%s)" % (domainObj.domain_id, domainObj.domain_description, "{:.1%}".format(domainObj.domain_protein_cov)))
            else:
                domain_string.append("None")
            outline.append("%s" % ";".join(domain_string))
            output.append("%s" % ",".join([str(field) for field in outline]))
        with open(matrix_f, 'w') as matrix_fh:
            matrix_fh.write("\n".join(output) + "\n")


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    target_cluster_id_f = args['--cluster_ids']
    outprefix = args['--outprefix']
    taxon_order_f = args['--taxon_order']
    cluster_metrics_f = args['--cluster_metrics']
    clustering_f = args['--clustering']
    DOMAIN_PROTEOME_COV = float(args['--domain_proteome_cov'])

    print "[+] Start ..."
    dataCollection = DataCollection(outprefix)
    dataCollection.parse_taxon_order_f(taxon_order_f)
    dataCollection.parse_target_cluster_id_f(target_cluster_id_f)
    dataCollection.parse_clustering_f(clustering_f)
    dataCollection.parse_cluster_metrics(cluster_metrics_f, DOMAIN_PROTEOME_COV)
    dataCollection.write_output()
