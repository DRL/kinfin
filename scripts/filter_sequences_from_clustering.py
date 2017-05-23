#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: filter_sequences_from_clustering.py      -c <FILE> [-o <STR>] [-i <FILE>] [-e <FILE>]
                                                    [-v] [-h|--help]

    Options:
        -h --help                           show this
        -o, --out_prefix <STR>              Output prefix for filtered BLAST files
        -c, --cluster_f <FILE>              Orthogroups
        -i, --include_id_f <FILE>           File containing headers of sequences to be included
        -e, --exclude_id_f <FILE>           File containing headers of sequences to be included
        -v, --verbose                       Verbose output
"""

from __future__ import division
from docopt import docopt
import os
import sys
from collections import defaultdict
from collections import Counter

class DataCollection():
    def __init__(self, include_ids, exclude_ids):
        self.include_ids = include_ids
        self.exclude_ids = exclude_ids
        self.status_by_protein_id = {}

    def set_status_by_protein_id(self):
        if self.exclude_ids:
            self.status_by_protein_id = defaultdict(lambda: True)
            for exclude_id in list(self.exclude_ids):
                self.status_by_protein_id[exclude_id] = False
        elif self.include_ids:
            self.status_by_protein_id = defaultdict(lambda: False)
            for include_id in list(self.include_ids):
                self.status_by_protein_id[include_id] = True
        else:
            pass

    def filter_cluster_f(self, cluster_f, out_prefix):
        print "[+] Filtering %s ..." % (cluster_f)
        included_out_f = "%s.included.txt" % (os.path.basename(cluster_f))
        excluded_out_f = "%s.excluded.txt" % (os.path.basename(cluster_f))
        statsfile = "%s.filtered.stats.txt" % (os.path.basename(cluster_f))
        if out_prefix:
            included_out_f = "%s.%s" % (out_prefix, included_out_f)
            excluded_out_f = "%s.%s" % (out_prefix, excluded_out_f)
            statsfile = "%s.%s" % (out_prefix, statsfile)
        cluster_stats_dict = {"total": 0, "excluded_singletons": 0, "excluded_non_singletons": 0, "included_singletons": 0, "included_non_singletons": 0}
        protein_stats_dict = {"total": 0, "excluded": 0, "included": 0}
        included_clusters = []
        excluded_clusters = []
        for line in read_file(cluster_f):
            cluster_id, protein_string = line.rstrip("\n").split(": ")
            cluster_stats_dict["total"] += 1
            protein_ids = protein_string.split(" ")
            protein_stats_dict["total"] += len(protein_ids)
            cluster_status = [self.status_by_protein_id[protein_id] for protein_id in protein_ids]
            singleton = True if len(protein_ids) == 1 else False
            if all(cluster_status):
                included_clusters.append(line)
                protein_stats_dict["included"] += len(protein_ids)
                if singleton is True:
                    cluster_stats_dict["included_singletons"] += 1
                else:
                    cluster_stats_dict["included_non_singletons"] += 1
            else:
                excluded_clusters.append(line)
                protein_stats_dict["excluded"] += len(protein_ids)
                if singleton is True:
                    cluster_stats_dict["excluded_singletons"] += 1
                else:
                    cluster_stats_dict["excluded_non_singletons"] += 1
        with open(included_out_f, 'w') as included_out_fh:
            included_out_fh.write("".join(included_clusters))
        with open(excluded_out_f, 'w') as excluded_out_fh:
            excluded_out_fh.write("".join(excluded_clusters))
        with open(statsfile, 'w') as stats_fh:
            stats_fh.write("file=%s;stats=clusters;total=%s;included_singletons=%s;included_non_singletons=%s;excluded_singletons=%s;excluded_non_singletons=%s\n" % (
                cluster_f,
                cluster_stats_dict["total"],
                cluster_stats_dict["included_singletons"],
                cluster_stats_dict["included_non_singletons"],
                cluster_stats_dict["excluded_singletons"],
                cluster_stats_dict["excluded_non_singletons"]))
            stats_fh.write("file=%s;stats=proteins;total=%s;included=%s;excluded=%s\n" % (
                cluster_f,
                protein_stats_dict["total"],
                protein_stats_dict["included"],
                protein_stats_dict["excluded"]))


def read_file(infile):
    if not os.path.exists(infile):
        sys.exit("[X] File %s does not exist." % (infile))
    else:
        with open(infile) as fh:
            for line in fh:
                yield line


def get_ids(include_id_f, exclude_id_f):
    include_ids = []
    exclude_ids = []
    if include_id_f and exclude_id_f:
        sys.exit("[X] Please provide an argument for EITHER --include_id_f or --exclude_id_f. Can't be both.")
    elif include_id_f:
        include_ids = parse_ids_f(include_id_f)
    elif exclude_id_f:
        exclude_ids = parse_ids_f(exclude_id_f)
    else:
        sys.exit("[X] Please provide an argument for --include_id_f or --exclude_id_f.")
    return set(include_ids), set(exclude_ids)


def parse_ids_f(id_f):
    ids = []
    for line in read_file(id_f):
        ids.append(line.rstrip("\n"))
    return set(ids)


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    out_prefix = args['--out_prefix']
    cluster_f = args['--cluster_f']
    include_id_f = args['--include_id_f']
    exclude_id_f = args['--exclude_id_f']
    verbose_flag = args['--verbose']

    # include/exclude
    include_ids, exclude_ids = get_ids(include_id_f, exclude_id_f)
    # Init data
    dataCollection = DataCollection(include_ids, exclude_ids)
    # set status of protein_ids depending on include_ids/exclude_ids
    dataCollection.set_status_by_protein_id()
    # filte BLAST based on sequence_ids
    dataCollection.filter_cluster_f(cluster_f, out_prefix)
