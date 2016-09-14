#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: kinfin-d.py      [-s <FILE>] [-g <FILE>] [-c <FILE>]
                        [-d <FILE>] [-b <FILE>] [--fasta_dir <DIR>]
                        [--nodesdb <FILE>] [--delimiter <STRING>]
                        [-f <FLOAT>] [-n <INT>] [--min <INT>] [--max <INT>]
                        [-l <INT>] [-r <INT>]
                        [--fontsize <INT>] [--plotsize INT,INT]
                        [-o <PREFIX>] [-p <PLOTFORMAT>]
                        [-h|--help]

    Options:
        -h --help                           show this

        Input files
            -s, --species_file <FILE>           SpeciesIDs.txt used in OrthoFinder
            -g, --groups <FILE>                 OrthologousGroups.txt produced by OrthoFinder
            -c, --classification_file <FILE>    SpeciesClassification file
            -d, --domain_file <FILE>            InterProScan TSV file
            -b, --blast_file <FILE>             BLAST results of similarity search against NR/Uniref90
            --fasta_dir <DIR>               Directory of FASTA files
            --nodesdb <FILE>                    nodesdb file (sames as blobtools nodesDB file)

        General options
            -o, --outprefix <STR>               Output prefix
            --delimiter <STRING>                Delimiter between proteome prefix and protein name [default: "."]

        "Fuzzy"-Orthology-groups
            -f <FLOAT>                          "Fuzzyness"-factor F [default: 0.9]
            -n <INT>                            Count of proteins in (F) of cluster [default: 1]
            --min <INT>                         Minimum count of proteins in (1-F) of cluster [default: 0]
            --max <INT>                         Maximum count of proteins in (1-F) of cluster[default: 100]

        Plotting
            -r, --repetitions <INT>             Number of repetitions for rarefaction curves [default: 30]
            --fontsize <INT>                    Fontsize for plots [default: 16]
            --plotsize <INT,INT>                Size (WIDTH,HEIGHT) for plots [default: 24,12]
            -p, --plotfmt <STR>                 Plot formats [default: png]

"""

from __future__ import division
import sys
from os.path import basename, isfile, abspath, splitext, join, exists
from os import getcwd, mkdir
from docopt import docopt, DocoptExit
from collections import defaultdict
from collections import Counter
from itertools import chain, combinations

########################################################################
# General constructs
########################################################################

tree_dict = lambda : defaultdict(tree_dict)

########################################################################
# General functions
########################################################################

def progress(iteration, steps, max_value):
    if int(iteration) == int(max_value):
        sys.stdout.write('\r')
        print "[PROGRESS]\t: %d%%" % (100)
    elif int(iteration) % int(steps+1) == 0:
        sys.stdout.write('\r')
        print "[PROGRESS]\t: %d%%" % (float(int(iteration)/int(max_value))*100),
        sys.stdout.flush()
    else:
        pass

def read_file(f):
    if not exists(f):
        sys.exit("[ERROR] - File %s does not exist." % (f))
    with open(f) as fh:
        for line in fh:
            yield line.rstrip("\n")

########################################################################
# Classes
########################################################################

class DomainRecordObj(object):
    """docstring for DomainObj"""
    def __init__(self, domain_id, protein_id, domain_source, domain_evalue, domain_desc, domain_start, domain_stop):
        self.id = domain_id
        self.protein_id = protein_id
        self.source = domain_source
        self.evalue = domain_evalue
        self.desc = domain_desc
        self.start = domain_start
        self.stop = domain_stop

class ProteinCollection():

class DomainCollection():
    def __init__(self):
        self.f = None
        self.by_protein_by_source = defaultdict(lambda: defaultdict(list))
        self.by_source_by_protein = defaultdict(lambda: defaultdict(list))
        self.by_source = {}
        self.counts = tree_dict()
        self.proteins_seen = set()

    def parse_domains(self, domain_f):
        '''
        parse Domains from InterProScan TSV output
        '''
        for line in read_file(domain_f):
            temp = line.split()
            protein_id = temp[0]
            domain_source = temp[3] # CDD, PIRSF, Pfam, Phobius, ProSiteProfiles, SMART, SUPERFAMILY, SignalP_EUK, TIGRFAM, TMHMM
            domain_id = temp[4]
            domain_evalue = temp[-3]
            domain_stop = temp[-4]
            domain_start = temp[-5]
            domain_desc = None
            if len(" ".join(temp[5:-5])):
                domain_desc = "\"%s\"" % " ".join(temp[5:-5])
            domain = DomainRecordObj(domain_id, protein_id, domain_source, domain_evalue, domain_desc, domain_start, domain_stop)
            self.add_domain(domain)
        self.f = domain_f

    def add_domain(self, domain):
        self.by_protein_by_source[domain.protein_id][domain.source].append(domain)
        self.by_source_by_protein[domain.source][domain.protein_id].append(domain)
        self.counts[domain.source][domain.id] = self.counts.get(domain.id, 0) + 1






if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    try:
        domain_f = args['--domain_file']
    except docopt.DocoptExit:
        print __doc__.strip()

    domainCollection = DomainCollection()
    print "[STATUS] - Parsing domains from %s" % (domain_f)
    domainCollection.parse_domains(domain_f)
    domainCollection.count_combinations()
