#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage:
    functional_annotation_of_clusters.py      all -f <FILE> -c <FILE> [-p <FLOAT>] [-x <FLOAT>] [-o <STRING>]
    functional_annotation_of_clusters.py      synapo -f <FILE> -c <FILE> -t <FILE> [-a <FILE>] [-n <FLOAT>] [-p <FLOAT>] [-x <FLOAT>] [-o <STRING>]
                                                        [-h|--help]

    Options:
        -h, --help                              show this
        -f, --cluster_domain_annotation <FILE>  cluster_domain_annotation.*.txt file (* = IPR/Pfam/GO/SignalP_Euk)
        -c, --cluster_counts_by_taxon <FILE>    Needed for getting number of proteins for clusters without domains
        -p, --domain_protein_cov <FLOAT>        Minimum protein coverage of domain in cluster [default: 0.75]
        -x, --domain_taxon_cov <FLOAT>          Minimum taxon coverage by proteins with domain in cluster [default: 0.75]

        -t, --tree_cluster_metrics <FILE>       tree.cluster_metrics.txt file
        -a, --node_names <FILE>                 CSV file of nodes and names (only output for those nodes will be supplied)
        -n, --node_taxon_cov <FLOAT>            Minimum taxon coverage of cluster in tree.cluster_metrics.txt [default: 0.75]

        -o, --outprefix <STRING>                Output prefix

"""
import sys
from docopt import docopt
import os


def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[ERROR] - File '%s' does not exist." % (infile))
    print("[+] Parsing %s ..." % (infile))
    with open(infile) as fh:
        for line in fh:
            line = line.replace(r'\r', '\n')
            if not line.startswith("#"):
                yield line.rstrip("\n")


def write_file(out_f, outprefix, header, strings):
    if outprefix:
        if outprefix.endswith("/"):
            if not os.path.exists(outprefix):
                os.mkdir(outprefix)
            out_f = "%s" % os.path.join(outprefix, out_f)
        else:
            out_f = "%s.%s" % (outprefix, out_f)
    print("[+] \t Writing file %s ..." % (out_f))
    with open(out_f, 'w') as out_fh:
        out_fh.write("%s\n" % (header))
        out_fh.write("%s\n" % "\n".join(strings))


class DomainObj():
    def __init__(self, domain_id, domain_description, domain_proteome_cov, domain_protein_cov):
        self.domain_id = domain_id
        self.domain_description = domain_description
        self.domain_proteome_cov = domain_proteome_cov
        self.domain_protein_cov = domain_protein_cov


class ClusterObj():
    def __init__(self, cluster_id, protein_count, taxon_count, boolean):
        self.cluster_id = cluster_id
        self.protein_count = protein_count
        self.taxon_count = taxon_count
        self.node_taxon_cov = 0.0
        self.status = boolean
        self.domainObjs = []

    def is_valid(self):
        return self.status

    def add_domainObj(self, domainObj):
        self.domainObjs.append(domainObj)

    def get_domain_line_all(self):
        line_cluster_domains = []
        line_cluster_domains.append(self.cluster_id)
        line_cluster_domains.append(self.protein_count)
        line_cluster_domains.append(self.taxon_count)
        if self.domainObjs:
            domain_ids = []
            domain_descriptions = []
            for domainObj in sorted(self.domainObjs, key=lambda x: x.domain_proteome_cov, reverse=True):
                domain_ids.append(domainObj.domain_id)
                domain_descriptions.append(domainObj.domain_description)
            line_cluster_domains.append(";".join(domain_ids))
            line_cluster_domains.append(";".join(domain_descriptions))
        else:
            line_cluster_domains.append("None")
            line_cluster_domains.append("None")
        return "\t".join([str(x) for x in line_cluster_domains])

    def get_domain_line_synapo(self, node_id, node_description):
        line_cluster_domains = []
        line_cluster_domains.append(self.cluster_id)
        line_cluster_domains.append(node_id)
        line_cluster_domains.append(node_description)
        line_cluster_domains.append(self.protein_count)
        line_cluster_domains.append(self.taxon_count)
        line_cluster_domains.append(self.node_taxon_coverage)
        if self.domainObjs:
            domain_ids = []
            domain_descriptions = []
            for domainObj in sorted(self.domainObjs, key=lambda x: x.domain_proteome_cov, reverse=True):
                domain_ids.append(domainObj.domain_id)
                domain_descriptions.append(domainObj.domain_description)
            line_cluster_domains.append(";".join(domain_ids))
            line_cluster_domains.append(";".join(domain_descriptions))
        else:
            line_cluster_domains.append("None")
            line_cluster_domains.append("None")
        return "\t".join([str(x) for x in line_cluster_domains])


class DataCollection():
    def __init__(self, args):
        self.all_flag = args['all']
        self.synapo_flag = args['synapo']
        self.cluster_metrics_f = args['--cluster_domain_annotation']
        self.cluster_counts_by_taxon_f = args['--cluster_counts_by_taxon']
        self.tree_cluster_metrics_f = args['--tree_cluster_metrics']
        self.node_names_f = args['--node_names']
        self.outprefix = args['--outprefix']
        self.NODE_TAXON_COV = self.get_float(args, '--node_taxon_cov')
        self.DOMAIN_PROTEIN_COV = self.get_float(args, '--domain_protein_cov')
        self.DOMAIN_TAXON_COV = self.get_float(args, '--domain_taxon_cov')
        self.cluster_ids = []
        self.clusterObjs_by_cluster_id = {}
        self.cluster_ids_by_node_id = {}
        self.node_ids = []
        self.node_description_by_node_id = {}
        if self.all_flag:
            self.parse_cluster_counts(True)  # all will be written
        elif self.synapo_flag:
            self.parse_cluster_counts(False)  # will be set to true for those that should be written
            self.parse_node_names()
            self.parse_tree_cluster_metrics()
        else:
            pass
        self.parse_cluster_metrics()
        self.write()

    def add_cluster_id_to_node_id(self, cluster_id, node_id):
        if node_id not in self.cluster_ids_by_node_id:
            self.node_ids.append(node_id)
            self.cluster_ids_by_node_id[node_id] = []
        self.cluster_ids_by_node_id[node_id].append(cluster_id)

    def parse_tree_cluster_metrics(self):
        for line in read_file(self.tree_cluster_metrics_f):
            col = line.split('\t')
            cluster_id = col[0]
            node_id = col[1]
            node_taxon_cov = float(col[3])
            if node_taxon_cov >= self.NODE_TAXON_COV:
                if self.node_description_by_node_id:
                    if node_id in self.node_description_by_node_id:
                        self.clusterObjs_by_cluster_id[cluster_id].status = True
                        self.clusterObjs_by_cluster_id[cluster_id].node_taxon_coverage = node_taxon_cov
                        self.add_cluster_id_to_node_id(cluster_id, node_id)
                else:
                    self.clusterObjs_by_cluster_id[cluster_id].status = True
                    self.clusterObjs_by_cluster_id[cluster_id].node_taxon_coverage = node_taxon_cov
                    self.add_cluster_id_to_node_id(cluster_id, node_id)

    def parse_node_names(self):
        if self.node_names_f:
            for line in read_file(self.node_names_f):
                col = line.split(',')
                try:
                    node_id = col[0]
                    node_description = col[1]
                except IndexError:
                    sys.exit("[ERROR] - File '%s' is not CSV or does not contain two columns." % (self.node_names_f))
                self.node_description_by_node_id[node_id] = node_description

    def parse_cluster_metrics(self):
        for line in read_file(self.cluster_metrics_f):
            col = line.split("\t")
            cluster_id = col[0]
            domain_id = col[2]
            domain_description = col[3]
            protein_count = int(col[4])
            protein_count_with_domain = int(col[5])
            domain_protein_cov = protein_count_with_domain / protein_count
            domain_proteome_cov = float(col[6])
            # print cluster_id, domain_id, protein_count, protein_count_with_domain, domain_protein_cov, domain_protein_cov
            if self.clusterObjs_by_cluster_id[cluster_id].is_valid():
                if domain_proteome_cov >= self.DOMAIN_TAXON_COV:
                    # print cluster_id, "PASS proteome_cov %s >= %s" % (domain_proteome_cov, self.DOMAIN_TAXON_COV)
                    if domain_protein_cov >= self.DOMAIN_PROTEIN_COV:
                        # print cluster_id, "PASS protein_cov %s >= %s" % (domain_protein_cov, self.DOMAIN_PROTEIN_COV)
                        domainObj = DomainObj(domain_id, domain_description, domain_proteome_cov, domain_protein_cov)
                        self.add_domainObj(cluster_id, domainObj)
                    # else:
                        # print cluster_id, "FAIL protein_cov %s < %s" % (domain_protein_cov, self.DOMAIN_PROTEIN_COV)
                # else:
                    # print cluster_id, "FAIL proteome_cov %s < %s" % (domain_proteome_cov, self.DOMAIN_TAXON_COV)

    def parse_cluster_counts(self, boolean):
        for line in read_file(self.cluster_counts_by_taxon_f):
            col = line.split()
            cluster_id = col[0]
            protein_count = sum([int(count) for count in col[1:]])
            taxon_count = sum([1 for count in col[1:] if int(count) > 0])
            self.clusterObjs_by_cluster_id[cluster_id] = ClusterObj(cluster_id, protein_count, taxon_count, boolean)
            self.cluster_ids.append(cluster_id)

    def get_float(self, args, key):
        try:
            value = float(args[key])
        except ValueError:
            sys.exit("[X] '%s' must be float" % (key))
        if value > 1 or value < 0:
            sys.exit("[X] '%s' must be between 0 and 1" % (key))
        return value

    def add_domainObj(self, cluster_id, domainObj):
        self.clusterObjs_by_cluster_id[cluster_id].add_domainObj(domainObj)

    def write(self):
        output = []
        header = []
        out_f = ''
        if self.all_flag:
            out_f = 'cluster_functional_annotation.all.p%s.x%s.tsv' % ("{0:.0f}".format(self.DOMAIN_PROTEIN_COV * 100), "{0:.0f}".format(self.DOMAIN_TAXON_COV * 100))
            header = "\t".join(["#cluster_id", "protein_count", "taxon_count", "domain_ids", "domain_description"])
            for cluster_id in self.cluster_ids:
                clusterObj = self.clusterObjs_by_cluster_id[cluster_id]
                line = clusterObj.get_domain_line_all()
                output.append(line)
        elif self.synapo_flag:
            out_f = 'cluster_functional_annotation.synapo.p%s.x%s.n%s.tsv' % ("{0:.0f}".format(self.DOMAIN_PROTEIN_COV * 100), "{0:.0f}".format(self.DOMAIN_TAXON_COV * 100), "{0:.0f}".format(self.NODE_TAXON_COV * 100))
            header = "\t".join(["#cluster_id", "node_id", "node_description", "protein_count", "taxon_count", "node_taxon_coverage", "domain_ids", "domain_description"])
            for node_id in self.node_ids:
                node_description = self.node_description_by_node_id.get(node_id, node_id)
                clusterObjs = [self.clusterObjs_by_cluster_id[cluster_id] for cluster_id in self.cluster_ids_by_node_id[node_id]]
                for clusterObj in sorted(clusterObjs, key=lambda x: x.node_taxon_cov, reverse=True):
                    line = clusterObj.get_domain_line_synapo(node_id, node_description)
                    output.append(line)
        write_file(out_f, self.outprefix, header, output)


if __name__ == "__main__":
    __version__ = 0.3
    args = docopt(__doc__)
    print("[+] Start ...")
    dataCollection = DataCollection(args)
