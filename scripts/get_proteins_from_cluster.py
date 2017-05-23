#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: get_proteins_from_cluster.py     -g <FILE> [--header <FILE>]
                                        [-c <STRING>] [--clusters <FILE>] [-s]
                                        [-o <STR>]
                                        [-h|--help]

    Options:
        -h --help                       show this
        -g, --groups <FILE>             OrthologousGroups.txt produced by OrthoFinder
        --header <FILE>                 Filter based on sequence IDs in file
        -c, --cluster <STRING>          Filter based on cluster ID
        --clusters <FILE>               Filter based on cluster IDs in file
        -o, --outprefix <STR>           Outprefix
        -s, --single_out_file           Write all proteins to a single file

"""

from __future__ import division
from docopt import docopt
import sys
from os.path import basename, isfile, abspath, splitext, join, exists

def parse_headers(header_f):
    headers = {}
    with open(header_f) as header_fh:
        for line in header_fh:
            header = line.rstrip("\n")
            if header in headers:
                sys.exit("[-] header %s repeated" % (header))
            else:
                headers[header] = None
    return headers

def parse_clusters(clusters):
    clusters = {}
    with open(cluster_f) as cluster_fh:
        for line in cluster_fh:
            cluster = line.rstrip("\n")
            if cluster in clusters:
                sys.exit("[-] cluster %s repeated" % (cluster))
            else:
                clusters[cluster] = None
    return clusters

def parse_groups(group_f):
    output = {}
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

def write_output(output, outprefix):
    headers_found = set([k for k, v in headers.iteritems() if v])
    clusters_found = set([k for k, v in clusters.iteritems() if v])
    if headers:
        print "[+] Found %s of headers ..." % "{:.0%}".format(len(headers_found)/len(headers))
    if clusters:
        print "[+] Found %s of clusters ..." % "{:.0%}".format(len(clusters_found)/len(clusters))
    stats_f = "%s.parse_stats.txt" % (splitext(basename(groups_f))[0])
    if outprefix:
        stats_f = "%s.%s" % (outprefix, stats_f)
    if headers_found or clusters_found:
        print "[+] Writing files ..."
        if not single_out_file:
            stats_lines = []
            for clusterID, proteins in output.items():
                protein_lines = []
                protein_lines += proteins
                proteins_total = len(proteins)
                proteins_target = len([x for x in proteins if x in headers_found])
                proteins_non_target = proteins_total - proteins_target
                out_f = "%s.%s.txt" % (splitext(basename(groups_f))[0], clusterID)
                if outprefix:
                    out_f = "%s.%s" % (outprefix, out_f)
                with open(out_f, 'w') as out_fh:
                    out_fh.write("\n".join(protein_lines))
                stats_lines.append("%s total=%s target=%s non-target=%s" % (clusterID, proteins_total, proteins_target, proteins_non_target))
            with open(stats_f, 'w') as stats_fh:
                stats_fh.write("\n".join(stats_lines))
        else:
            protein_lines = []
            out_f = "%s.protein_ids.txt" % (splitext(basename(cluster_f))[0])
            if outprefix:
                out_f = "%s.%s" % (outprefix, out_f)
            for clusterID, proteins in output.items():
                protein_lines += proteins
            with open(out_f, 'w') as out_fh:
                out_fh.write("\n".join(protein_lines))


if __name__ == "__main__":
    __version__ = 0.2
    args = docopt(__doc__)
    groups_f = args['--groups']
    header_f = args['--header']
    cluster_id = args['--cluster']
    cluster_f = args['--clusters']
    single_out_file = args['--single_out_file']
    outprefix = args['--outprefix']

    headers = {}
    clusters = {}

    print "[+] Start ..."
    if header_f:
        print "[+] Parsing headers in %s ..." % header_f
        parse_type = 'header'
        headers = parse_headers(header_f)
    elif cluster_id:
        print "[+] Getting cluster %s ..." % cluster_id
        clusters[cluster_id] = None
    elif cluster_f:
        print "[+] Parsing clusters in %s ..." % cluster_f
        parse_type = 'cluster'
        clusters = parse_clusters(cluster_f)
    else:
        sys.exit(__doc__.strip())
    print "[+] Parse groups %s ..." % groups_f
    output = parse_groups(groups_f)
    write_output(output, outprefix)
