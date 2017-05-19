#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: compare_clusterings.py                 -c <FILE> -d <DIR> -s <FILE> -t <FILE> [-o <STRING>]
                                              [-h|--help]

    Options:
        -h --help                           show this
        -c, --config <FILE>                 Config file
        -d, --clustering_dir <DIR>          Clustering dir
        -s, --sequence <FILE>               SequenceID.txt file
        -t, --taxon <FILE>                  Taxon file
        -o, --outprefix <STRING>            Output prefix [default: results]

"""

import re
import sys
import operator
from docopt import docopt
from collections import deque
import os


#class BarChartObj():
#    def __init__(self, IV):
#        self.IV = IV
#        self.proteins_conserved = {'non-singletons-split' : set(), 'non-singletons' : set(), 'singletons' : set()}
#        self.proteins_non_conserved = {'non-singletons-split' : set(), 'non-singletons' : set(), 'singletons' : set()}
#        self.proteins_singletons = {'non-singletons-split' : set(), 'non-singletons' : set(), 'singletons' : set()}
#
#
#class DataCollection():
#    def __init__(self):
#        self.clusterCollections = []
#        self.protein_ids_by_hierarchality_pattern = {} # tuples of hierarchality
#        self.hierarchality_pattern_by_protein_id = {}
#        self.hierarchality_pattern_by_protein_ids = {}
#        self.protein_ids_trajectory = {}
#        self.cluster_ids_trajectory = {}
#
#
#    def add_to_hierarchy(self, protein_ids_hierarchy, protein_ids):
#        if protein_ids_hierarchy[protein_ids]:
#            protein_ids_hierarchy[protein_ids]
#
#    def add_clusterCollection(self, clusterCollection):
#        barChartObj = BarChartObj(clusterCollection.IV)
#
#        cluster_counts_by_split = {'conserved': {'singleton': 0, 'non_singleton' : 0}, 'hierarchical' : {'singleton': 0, 'non_singleton' : 0}, 'non_hierarchical' : 0}
#        protein_counts_by_split = {'conserved': {'singleton': 0, 'non_singleton' : 0}, 'hierarchical' : {'singleton': 0, 'non_singleton' : 0}, 'non_hierarchical' : 0}
#        cluster_ids_by_split = {'conserved': {'singleton': set(), 'non_singleton' : set()}, 'hierarchical' : {'singleton': set(), 'non_singleton' : set()}, 'non_hierarchical' : set()}
#        protein_ids_by_split = {'conserved': {'singleton': set(), 'non_singleton' : set()}, 'hierarchical' : {'singleton': set(), 'non_singleton' : set()}, 'non_hierarchical' : set()}
#        if self.clusterCollections == []:
#            self.clusterCollections.append(clusterCollection)
#            cluster_counts_by_split['conserved']['singleton'] += clusterCollection.cluster_singleton_count
#            cluster_counts_by_split['conserved']['non_singleton'] += clusterCollection.cluster_non_singleton_count
#            protein_counts_by_split['conserved']['singleton'] += clusterCollection.protein_singleton_count
#            protein_counts_by_split['conserved']['non_singleton'] += clusterCollection.protein_count
#        else:
#            parent_clusterCollection = self.clusterCollections[-1]
#            child_clusterCollection = child_clusterCollection
#            non_hierarchical_cluster_ids = set()
#
#            for parent_protein_ids, prev_cluster_id in parent_clusterCollection.cluster_id_by_protein_ids.items():
#                # for each previous cluster
#                if parent_protein_ids in child_clusterCollection.cluster_id_by_protein_ids:
#                    # if cluster is present in current clustering => conserved, no-split
#                    split = 'None'
#                else:
#                    # cluster is split somehow: hierarchical (singleton/non-singleton), non-hierarchical (non-singleton)
#                    protein_ids_split_in_singleton = []
#                    protein_ids_split_in_non_singleton = []
#                    for parent_protein_id in parent_protein_ids:
#                        # for each protein in previous cluster
#                        if parent_protein_id in child_clusterCollection.cluster_id_by_protein_ids:
#                            # if singleton
#                            protein_ids_split_in_singleton.append(parent_protein_id)
#                        else:
#                            # not a singleton
#                            protein_ids_split_in_non_singleton.append(parent_protein_id)
#                        if protein_ids_split_in_non_singleton:
#                            non_singleton_clusters = set()
#                            for protein_id in protein_ids_split_in_non_singleton:
#                                non_singleton_clusters.add(clusterCollection.cluster_id_by_protein_id[protein_id])
#                            for non_singleton_cluster in non_singleton_clusters:
#                                # check if additional proteins are present
#                                protein_ids_in_non_singleton_cluster = clusterCollection.clusterObj_by_cluster_id[non_singleton_cluster]
#                                additional_protein_ids = protein_ids_in_non_singleton_cluster.difference(parent_protein_ids)
#                                if additional_protein_ids:
#                                    if not non_singleton_cluster in non_hierarchical_cluster_ids:
#                                        non_hierarchical_cluster_count += 1
#                                        non_hierarchical_cluster_ids.add(non_singleton_cluster)
#                                    non_hierarchical_split_count += 1
#                                else:
#                                    hierarchical_cluster_count += 1
#                                    hierarchical_split_count += 1
#                                    "non-hierarchical"
#                else:
#                    if len(protein_ids) == 1:
#                        hierarchality = "split_singleton"
#                    else:
#                        hierarchality = "split_non_singleton"
#                for protein_id in protein_ids:
#                    if not protein_id in self.hierarchality_pattern_by_protein_id:
#                        self.hierarchality_pattern_by_protein_id[protein_id] = []
#                    self.hierarchality_pattern_by_protein_id[protein_id].append(hierarchality)
#                if not protein_ids in self.hierarchality_pattern_by_protein_ids:
#                    self.hierarchality_pattern_by_protein_ids[protein_ids] = []
#                self.hierarchality_pattern_by_protein_ids[protein_ids].append(hierarchality)
#            self.clusterCollections.append(clusterCollection)
#
#    def generate_output(self):
#        for protein_ids, hierarchality_patterns in self.hierarchality_pattern_by_protein_ids.items():
#            print "proteins: %s" % (len(protein_ids)), hierarchality_patterns
#


##########################################################################################################################


def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        for line in fh:
            yield line.rstrip("\n")


class CountObj():
    def __init__(self):
        self.singleton_clusters = 0
        self.singleton_proteins = 0
        self.monotaxon_clusters = 0
        self.monotaxon_proteins = 0
        self.multitaxon_clusters = 0
        self.multitaxon_proteins = 0
        self.clusters = 0
        self.proteins = 0

    def count_clusterObj(self, clusterObj):
        self.clusters += 1
        self.proteins += clusterObj.protein_count
        if clusterObj.cluster_type == 'singleton':
            self.singleton_clusters += 1
            self.singleton_proteins += 1
        elif clusterObj.cluster_type == 'monotaxon':
            self.monotaxon_clusters += 1
            self.monotaxon_proteins += clusterObj.protein_count
        elif clusterObj.cluster_type == 'multitaxon':
            self.multitaxon_clusters += 1
            self.multitaxon_proteins += clusterObj.protein_count
        else:
            sys.exit("[-] cluster %s has cluster_type %s" % clusterObj.cluster_id, clusterObj.cluster_type)


class ClusterObj():
    def __init__(self, IV, cluster_id, protein_ids):
        self.IV = IV
        self.cluster_id = cluster_id
        self.protein_ids = protein_ids
        self.protein_count = len(protein_ids)
        self.taxon_count = self.get_taxon_count(protein_ids)
        self.cluster_type = self.get_cluster_type()

    def get_taxon_count(self, protein_ids):
        taxon_list = []
        for protein_id in protein_ids:
            taxon_id = dataCollection.taxon_id_by_protein_id[protein_id]
            taxon_list.append(taxon_id)
        return len(set(taxon_list))

    def get_cluster_type(self):
        if self.taxon_count > 1:
            return "multitaxon"
        elif self.protein_count > 1:
            return "monotaxon"
        elif self.protein_count == 1:
            return "singleton"
        else:
            return None


class ClusterCollection():  # This is equivalent to orthogroup
    def __init__(self, clustering_id, clustering_f):
        self.clustering_id = clustering_id
        self.clustering_f = clustering_f
        self.countObj = CountObj()
        self.cluster_ids = []
        self.cluster_id_by_protein_id = {}
        self.cluster_id_by_protein_ids = {}
        self.clusterObj_by_cluster_id = {}

    def parse_clustering_f(self):
        print "[+] \t Parsing cluster file %s ..." % (self.clustering_f)
        for line in read_file(self.clustering_f):
            cluster_id, protein_string = line.rstrip("\n").split(": ")
            cluster_id = "%s|%s" % (self.clustering_id, cluster_id)        # add clustering_id as prefix to cluster name
            protein_ids = protein_string.split(" ")
            clusterObj = ClusterObj(self.clustering_id, cluster_id, protein_ids)
            self.countObj.count_clusterObj(clusterObj)
            for protein_id in clusterObj.protein_ids:
                self.cluster_id_by_protein_id[protein_id] = clusterObj.cluster_id
            self.cluster_id_by_protein_ids[frozenset(clusterObj.protein_ids)] = clusterObj.cluster_id
            self.cluster_ids.append(clusterObj.cluster_id)
            self.clusterObj_by_cluster_id[clusterObj.cluster_id] = clusterObj
            self.cluster_ids.append(clusterObj.cluster_id)


class DataCollection():
    def __init__(self, clustering_dir, out_prefix):
        self.clustering_dir = clustering_dir
        self.out_prefix = out_prefix

        self.taxon_id_by_taxon_idx = {}
        self.clustering_f_by_clustering_id = {}
        self.clustering_ids = []
        self.clusterCollection_by_clustering_id = {}
        self.taxon_id_by_protein_id = {}    # dict based on SequenceIDs.txt
        self.metaClusterObjs_by_metacluster_id = {}
        self.metacluster_id_by_cluster_id = {}
        self.metacluster_id_by_protein_id = {}

    def parse_taxon_f(self, taxon_f):
        print "[+] Parsing taxon file %s ..." % (taxon_f)
        for line in read_file(taxon_f):
            if not line.startswith("#"):
                try:
                    taxon_idx, taxon_id = line.split(",")
                except ValueError:
                    sys.exit("[-] Taxon file %s does not seem to be in the right format. Should have two fields separated by \",\"" % taxon_f)
                self.taxon_id_by_taxon_idx[taxon_id] = taxon_idx

    def parse_sequence_id_f(self, sequence_id_f):
        print "[+] Parsing SequenceID file %s ..." % (sequence_id_f)
        for line in read_file(sequence_id_f):
            try:
                idx_string, protein_string = line.split(": ")
                protein_id = protein_string.replace(":", "_").replace(",", "_")
                taxon_id = idx_string.split("_")[0]
            except ValueError:
                sys.exit("[-] SequenceID file %s does not seem to be in the right format. Should have two fields separated by \": \"" % sequence_id_f)
            self.taxon_id_by_protein_id[protein_id] = taxon_id

    def parse_config_f(self, config_f):
        print "[+] Parsing config file %s ..." % (config_f)
        for line in read_file(config_f):
            if not line.startswith("#"):
                try:
                    clustering_id, clustering_f = line.split(",")
                except ValueError:
                    sys.exit("[-] Config file %s does not seem to be in the right format. Should have two fields separated by \",\"" % config_f)
                self.clustering_ids.append(clustering_id)
                self.clustering_f_by_clustering_id[clustering_id] = clustering_f

    def generate_clusterCollections(self):
        print "[+] Generating clusterCollections ..."
        for clustering_id in self.clustering_ids:
            clustering_f = self.clustering_f_by_clustering_id[clustering_id]
            clustering_path = os.path.join(self.clustering_dir, clustering_f)
            clusterCollection = ClusterCollection(clustering_id, clustering_path)
            clusterCollection.parse_clustering_f()
            self.clusterCollection_by_clustering_id[clustering_id] = clusterCollection

    def get_cluster_ids_for_protein_id(self, protein_id):   # example how to access cluster_ids for protein_id
        cluster_ids = []
        for clustering_id in self.clustering_ids:
            clusterCollection = self.clusterCollection_by_clustering_id[clustering_id]
            cluster_id = clusterCollection.cluster_id_by_protein_id[protein_id]
            cluster_ids.append(cluster_id)
        return cluster_ids

    def get_protein_ids_for_cluster_id(self, cluster_id):
        clustering_id = cluster_id.split("|")[0]
        return self.clusterCollection_by_clustering_id[clustering_id].clusterObj_by_cluster_id[cluster_id].protein_ids

    def get_clusterObj_for_cluster_id(self, cluster_id):
        clustering_id = cluster_id.split("|")[0]
        return self.clusterCollection_by_clustering_id[clustering_id].clusterObj_by_cluster_id[cluster_id]

    def traverse(self, cluster_ids_seen, protein_ids_seen):
        print "# traverse"
        print "proteins seen", len(protein_ids_seen)
        print "clusters seen", len(cluster_ids_seen)
        if not all([value for key, value in cluster_ids_seen.items()]) or not all([value for key, value in protein_ids_seen.items()]):
            for protein_id, status in protein_ids_seen.items():
                if status is False:
                    protein_ids_seen[protein_id] = True
                    for cluster_id in self.get_cluster_ids_for_protein_id(protein_id):
                        if cluster_id not in cluster_ids_seen:
                            cluster_ids_seen[cluster_id] = False
            for cluster_id, status in cluster_ids_seen.items():
                if status is False:
                    cluster_ids_seen[cluster_id] = True
                    for protein_id in self.get_protein_ids_for_cluster_id(cluster_id):
                        if protein_id not in protein_ids_seen:
                            protein_ids_seen[protein_id] = False
            self.traverse(cluster_ids_seen, protein_ids_seen)
        else:
            print "Done"
            print len(cluster_ids_seen)
            print len(protein_ids_seen)
        return cluster_ids_seen, protein_ids_seen

    def generate_metaClusters(self):
        print "[+] Generating MetaClusters ..."
        idx = -1
        for clustering_id in self.clustering_ids:
            clusterCollection = self.clusterCollection_by_clustering_id[clustering_id]
            for cluster_id in clusterCollection.cluster_ids:
                if cluster_id not in self.metacluster_id_by_cluster_id:  # Novel cluster
                    idx += 1

                    metaClusterObj = MetaClusterObj(idx)
                    print "#", metaClusterObj.metacluster_id
                    cluster_ids_seen = {}
                    protein_ids_seen = {}
                    cluster_ids_seen[cluster_id] = True
                    for protein_id in self.get_protein_ids_for_cluster_id(cluster_id):
                        protein_ids_seen[protein_id] = False
                    cluster_ids_seen, protein_ids_seen = self.traverse(cluster_ids_seen, protein_ids_seen)
                    for cluster_id_seen in cluster_ids_seen:
                        self.metacluster_id_by_cluster_id[cluster_id_seen] = metaClusterObj.metacluster_id
        print self.metacluster_id_by_cluster_id
                    #    clusterObj = self.get_clusterObj_for_cluster_id(cluster_id_seen)
                    #    metaClusterObj.add_clusterObj(clusterObj)
                    #for protein_id_seen in protein_ids_seen:
                    #    self.metacluster_id_by_protein_id[protein_id_seen] = metaClusterObj.metacluster_id


class MetaClusterObj():
    def __init__(self, idx):
        self.metacluster_id = "MG%s" % idx
        self.cluster_ids_by_clustering_id = {}
        self.countObj_by_clustering_id = {}
        self.countObj = CountObj()

    def add_clusterObj(self, clusterObj):
        clustering_id = clusterObj.clustering_id
        if clustering_id not in self.cluster_ids_by_clustering_id:
            self.cluster_ids_by_clustering_id[clustering_id] = []
            self.countObj_by_clustering_id[clustering_id] = CountObj()
        self.cluster_ids_by_clustering_id[clustering_id].append(cluster_id)
        self.countObj_by_clustering_id[clustering_id].count_clusterObj(clusterObj)
        self.countObj.count_clusterObj(clusterObj)

# Questions
# how many metaclusters are there? what is the size distribution?
# how many clusters are stable across all/some inflations values?
#    - write cluster_ids that are equivalent
#    - What are their domain composition? GO terms?
#    - they must be sequence similarity islands
# are isoforms more stable than non-isoforms? => needs isoform input file
# what is the rate of splits (hierarchical/non-hierarchical)?
#   - how does that change across IVs?
#   - how many clusters split non-hierachically ?


if __name__ == "__main__":
    __version__ = 0.2
    args = docopt(__doc__)

    config_f = args['--config']
    taxon_f = args['--taxon']
    sequence_id_f = args['--sequence']
    clustering_dir = args['--clustering_dir']
    out_prefix = args['--outprefix']

    print "[+] Start ..."
    dataCollection = DataCollection(clustering_dir, out_prefix)
    dataCollection.parse_config_f(config_f)
    dataCollection.parse_taxon_f(taxon_f)
    dataCollection.parse_sequence_id_f(sequence_id_f)
    dataCollection.generate_clusterCollections()
    dataCollection.generate_metaClusters()
