#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: find_isoforms.py                      -i <FILE> [-o <STRING>]
                                              [-h|--help]

    Options:
        -h --help                           show this
        -i, --infile <FILE>                 Infile
        -o, --outprefix <STRING>            Output prefix

"""

import re
import sys
import operator
from docopt import docopt
from os.path import isfile, join, exists, realpath, dirname, basename

def parse_groups(group_f):
    with open(groups_f) as group_fh:
        for line in group_fh:
            clusterID, protein_string = line.rstrip("\n").split(": ")
            proteins = protein_string.split(" ")
            if headers:
                for protein in proteins:
                    if protein in headers:
                        if headers[protein] == None:
                            headers[protein] = clusterID
                            output[clusterID] = proteins
                        else:
                            sys.exit("[-] protein %s found more than once" % protein)
            else:
                if clusterID in clusters:
                    if clusters[clusterID] == None:
                        clusters[clusterID] = clusterID
                        output[clusterID] = proteins
                    else:
                        sys.exit("[-] cluster %s found more than once" % clusterID)
    return output

def read_file(infile):
    if not infile or not exists(infile):
        sys.exit("[ERROR] - File '%s' does not exist." % (infile))
    if infile.endswith(".gz"):
        import gzip
        with gzip.open(infile) as fh:
            for line in fh:
                yield line.rstrip("\n")
    else:
        with open(infile) as fh:
            for line in fh:
                yield line.rstrip("\n")

class BarChartObj():
    def __init__(self, IV):
        self.IV = IV
        self.proteins_conserved = {'non-singletons-split' : set(), 'non-singletons' : set(), 'singletons' : set()}
        self.proteins_non_conserved = {'non-singletons-split' : set(), 'non-singletons' : set(), 'singletons' : set()}
        self.proteins_singletons = {'non-singletons-split' : set(), 'non-singletons' : set(), 'singletons' : set()}


class DataCollection():
    def __init__(self):
        self.clusterCollections = []
        self.protein_ids_by_hierarchality_pattern = {} # tuples of hierarchality
        self.hierarchality_pattern_by_protein_id = {}
        self.hierarchality_pattern_by_protein_ids = {}
        self.protein_ids_trajectory = {}
        self.cluster_ids_trajectory = {}


    def add_to_hierarchy(self, protein_ids_hierarchy, protein_ids):
        if protein_ids_hierarchy[protein_ids]:
            protein_ids_hierarchy[protein_ids]

    def add_clusterCollection(self, clusterCollection):
        barChartObj = BarChartObj(clusterCollection.IV)

        cluster_counts_by_split = {'conserved': {'singleton': 0, 'non_singleton' : 0}, 'hierarchical' : {'singleton': 0, 'non_singleton' : 0}, 'non_hierarchical' : 0}
        protein_counts_by_split = {'conserved': {'singleton': 0, 'non_singleton' : 0}, 'hierarchical' : {'singleton': 0, 'non_singleton' : 0}, 'non_hierarchical' : 0}
        cluster_ids_by_split = {'conserved': {'singleton': set(), 'non_singleton' : set()}, 'hierarchical' : {'singleton': set(), 'non_singleton' : set()}, 'non_hierarchical' : set()}
        protein_ids_by_split = {'conserved': {'singleton': set(), 'non_singleton' : set()}, 'hierarchical' : {'singleton': set(), 'non_singleton' : set()}, 'non_hierarchical' : set()}
        if self.clusterCollections == []:
            self.clusterCollections.append(clusterCollection)
            cluster_counts_by_split['conserved']['singleton'] += clusterCollection.cluster_singleton_count
            cluster_counts_by_split['conserved']['non_singleton'] += clusterCollection.cluster_non_singleton_count
            protein_counts_by_split['conserved']['singleton'] += clusterCollection.protein_singleton_count
            protein_counts_by_split['conserved']['non_singleton'] += clusterCollection.protein_count
        else:
            parent_clusterCollection = self.clusterCollections[-1]
            child_clusterCollection = child_clusterCollection
            non_hierarchical_cluster_ids = set()

            for parent_protein_ids, prev_cluster_id in parent_clusterCollection.cluster_id_by_protein_ids.items():
                # for each previous cluster
                if parent_protein_ids in child_clusterCollection.cluster_id_by_protein_ids:
                    # if cluster is present in current clustering => conserved, no-split
                    split = 'None'
                else:
                    # cluster is split somehow: hierarchical (singleton/non-singleton), non-hierarchical (non-singleton)
                    protein_ids_split_in_singleton = []
                    protein_ids_split_in_non_singleton = []
                    for parent_protein_id in parent_protein_ids:
                        # for each protein in previous cluster
                        if parent_protein_id in child_clusterCollection.cluster_id_by_protein_ids:
                            # if singleton
                            protein_ids_split_in_singleton.append(parent_protein_id)
                        else:
                            # not a singleton
                            protein_ids_split_in_non_singleton.append(parent_protein_id)
                        if protein_ids_split_in_non_singleton:
                            non_singleton_clusters = set()
                            for protein_id in protein_ids_split_in_non_singleton:
                                non_singleton_clusters.add(clusterCollection.cluster_id_by_protein_id[protein_id])
                            for non_singleton_cluster in non_singleton_clusters:
                                # check if additional proteins are present
                                protein_ids_in_non_singleton_cluster = clusterCollection.clusterObjs_by_cluster_id[non_singleton_cluster]
                                additional_protein_ids = protein_ids_in_non_singleton_cluster.difference(parent_protein_ids)
                                if additional_protein_ids:
                                    if not non_singleton_cluster in non_hierarchical_cluster_ids:
                                        non_hierarchical_cluster_count += 1
                                        non_hierarchical_cluster_ids.add(non_singleton_cluster)
                                    non_hierarchical_split_count += 1
                                else:
                                    hierarchical_cluster_count += 1
                                    hierarchical_split_count += 1

                                    "non-hierarchical"




                else:
                    if len(protein_ids) == 1:
                        hierarchality = "split_singleton"
                    else:
                        hierarchality = "split_non_singleton"
                for protein_id in protein_ids:
                    if not protein_id in self.hierarchality_pattern_by_protein_id:
                        self.hierarchality_pattern_by_protein_id[protein_id] = []
                    self.hierarchality_pattern_by_protein_id[protein_id].append(hierarchality)
                if not protein_ids in self.hierarchality_pattern_by_protein_ids:
                    self.hierarchality_pattern_by_protein_ids[protein_ids] = []
                self.hierarchality_pattern_by_protein_ids[protein_ids].append(hierarchality)
            self.clusterCollections.append(clusterCollection)

    def generate_output(self):
        for protein_ids, hierarchality_patterns in self.hierarchality_pattern_by_protein_ids.items():
            print "proteins: %s" % (len(protein_ids)), hierarchality_patterns

class ClusterObj():
    def __init__(self, IV, cluster_id, protein_ids):
        self.IV = IV
        self.cluster_id = cluster_id
        self.protein_ids = protein_ids
        self.protein_count = len(protein_ids)
        self.singleton = True if len(protein_ids) == 1 else False
        self.validated = False

class ClusterCollection():
    def __init__(self, IV, cluster_f):
        self.IV = IV
        self.cluster_f = cluster_f

        self.cluster_id_by_protein_id = {}
        self.cluster_id_by_protein_ids = {}
        self.clusterObjs_by_cluster_id = {}

        self.cluster_ids_by_cluster_status = {}

        self.cluster_count = 0
        self.protein_count = 0
        self.cluster_singleton_count = 0
        self.protein_singleton_count = 0
        self.cluster_non_singleton_count = 0
        self.protein_non_singleton_count = 0

    def parse_cluster_f(self):
        print "[+] Parsing cluster file %s ..." % (self.cluster_f)
        for line in read_file(self.cluster_f):
            cluster_id_string, protein_string = line.rstrip("\n").split(": ")
            cluster_id = "%s:%s" % (self.IV, cluster_id_string)
            protein_ids = protein_string.split(" ")
            clusterObj = ClusterObj(self.IV, cluster_id, protein_ids)
            # count
            if clusterObj.singleton == True:
                self.cluster_singleton_count += 1
                self.protein_singleton_count += len(clusterObj.protein_count)
            else:
                self.cluster_non_singleton_count += 1
                self.protein_non_singleton_count += len(clusterObj.protein_count)
            self.cluster_count += 1
            self.protein_count += len(clusterObj.protein_count)
            # add
            for protein_id in clusterObj.protein_ids:
                self.cluster_id_by_protein_id[protein_id] = clusterObj.cluster_id
            self.cluster_id_by_protein_ids[frozenset(clusterObj.protein_ids)] = clusterObj.cluster_id
            self.clusterObjs_by_cluster_id[clusterObj.cluster_id] = clusterObj

class InputObj():
    def __init__(self, input_f):
        self.input_file_by_IV = self.parse_input_f(input_f)

    def parse_input_f(self, input_f):
        print "[+] Parsing input file %s ..." % (input_f)
        input_file_by_IV = {}
        for line in read_file(input_f):
            temp = line.split()
            IV = temp[0]
            cluster_f = temp[1]
            input_file_by_IV[IV] = cluster_f
        return input_file_by_IV

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    input_f = args['--infile']
    out_prefix = args['--outprefix']

    inputObj = InputObj(input_f)
    print "[+] Start ..."
    dataCollection = DataCollection()
    for IV, cluster_f in sorted(inputObj.input_file_by_IV.items()):
        clusterCollection = ClusterCollection(IV, cluster_f)
        clusterCollection.parse_cluster_f()
        dataCollection.add_clusterCollection(clusterCollection)
    dataCollection.generate_output()
