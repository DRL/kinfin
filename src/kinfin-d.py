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

def parse_species_ids(species_ids_f):
    proteome_fastas = {}
    for line in read_file(species_ids_f):
        number, proteome_fasta = line.split(": ")
        proteome_fastas[number] = proteome_fasta
    return proteome_fastas

def parse_nodesdb(nodesdb_f):
    nodesdb = {}
    nodesdb_count = 0
    nodes_count = 0
    for line in read_file(nodesdb_f):
        if line.startswith("#"):
            nodesdb_count = int(line.lstrip("# nodes_count = ").rstrip("\n"))
        else:
            nodes_count += 1
            node, rank, name, parent = line.rstrip("\n").split("\t")
            nodesdb[node] = {'rank' : rank, 'name' : name, 'parent' : parent}
            if (nodesdb_count):
                progress(nodes_count, 1000, nodesdb_count)
    return nodesdb

########################################################################
# Classes
########################################################################

class AloCollection():
    def __init__(self):
        self.proteomes = set() # set of proteome_idxs for orthofinder
        self.attributes = [] # list of attribute
        self.levels_by_attributes_by_proteome = {}

    def parse_classification(self, species_classification_f):
        '''
        species_classification_f : user defined CSV config file used for creating ALOs (Attribute-Level-Objects)
            - header: starts with '#', required fields: 'proteome_idx' (used for linking Orthofinder-Proteome/ProteinIDs)
        '''
        print "[STATUS] - Parsing SpeciesClassification file %s" % (species_classification_f)
        attributes = None
        levels_by_attributes_by_proteome = {}
        # reading species_classification_f
        for line in read_file(species_classification_f):
            if line.startswith("#"):
                attributes = [x.strip() for x in line.lstrip("#").split(",")]
                if not 'proteome_idx' in set(attributes):
                    sys.exit("[ERROR] - One of the attributes has to be 'proteome_idx', so that Orthofinder SequenceIDs and SpeciesIDs can be linked.\n\t%s" % (attributes))
            elif line.strip():
                temp = line.split(",")
                if not len(temp) == len(attributes):
                    sys.exit("[ERROR] - number of columns in line differs from header\n\t%s\n\t%s" % (attributes, temp))
                if temp[0] in self.proteomes:
                    sys.exit("[ERROR] - 'proteome_idx' should be unique. %s was encountered multiple times" % (temp[0]))
                proteome_idx = temp[0]
                self.proteomes.add(proteome)
                levels_by_attributes_by_proteome[proteome] = {x : '' for x in attributes}
                for idx, level in enumerate(temp):
                    attribute = attributes[idx]
                    levels_by_attributes_by_proteome[proteome][attribute] = level
            else:
                pass
        # Save
        self.attributes = attributes
        self.levels_by_attributes_by_proteome = levels_by_attributes_by_proteome

        # Add taxonomy if needed
        if 'taxid' in set(self.attributes):
            print "[+] - Attribute 'taxid' found, inferring taxonomic ranks from nodesDB..."
            self.update_classification()

        self.levels_by_attribute = {attribute : set() for attribute in self.attributes}
        self.proteome_by_level_by_attribute = {attribute : {} for attribute in self.attribute}

        for proteome in self.levels_by_attributes_by_proteome:
            for attribute in self.levels_by_attributes_by_proteome[proteome]:
                level = self.levels_by_attributes_by_proteome[proteome][attribute]
                if not level in self.proteome_by_level_by_attribute[attribute]:
                    self.proteome_by_level_by_attribute[attribute][level] = set()

        # new levels
        for proteomeID in levelIDs_by_rankID_by_proteomeID:
            for rankID in levelIDs_by_rankID_by_proteomeID[proteomeID]:
                levelID = levelIDs_by_rankID_by_proteomeID[proteomeID][rankID]
                if not levelID in self.proteomeIDs_by_levelID_by_rankID[rankID]:
                    self.proteomeIDs_by_levelID_by_rankID[rankID][levelID] = set()
                if not rankID in self.count_by_levelID_by_rankID:
                    self.count_by_levelID_by_rankID[rankID] = {}
                if not levelID in self.count_by_levelID_by_rankID[rankID]:
                    self.count_by_levelID_by_rankID[rankID][levelID] = 1
                else:
                    self.count_by_levelID_by_rankID[rankID][levelID] += 1
                self.levelIDs_by_rankID[rankID].add(levelID)
                self.levelIDs_by_rankID_by_proteomeID[proteomeID][rankID] = levelID
                self.proteomeIDs_by_levelID_by_rankID[rankID][levelID].add(proteomeID)
        self.proteomeIDs_count = len(self.proteomeIDs)
        self.levelIDs_count_by_rankID = {rankID : len(levels) for rankID, levels in self.levelIDs_by_rankID.items()}

    def update_classification(self):
        if not nodesdb_f:
            sys.exit("[ERROR] - Please specify a nodesDB file, so that taxonomic ranks can be inferred based on taxIDs\n")
        print "[+] - Parsing nodesDB %s" % (nodesdb_f)
        NODESDB = parse_nodesdb(nodesdb_f)
        for proteome in self.levels_by_attributes_by_proteome:
            taxid = self.levels_by_attributes_by_proteome[proteome]['taxid']
            lineage = get_lineage(taxid)
            # add lineage attribute/levels
            for taxrank in TAXRANKS:
                self.levels_by_attributes_by_proteome[proteome][taxrank] = lineage[taxrank]
            # remove taxid-levels
            del self.levels_by_attributes_by_proteome[proteome]['taxid']
        # remove taxid-attribute
        self.attributes.remove('taxid')
        # add taxranks to rank
        for taxrank in TAXRANKS:
            self.attributes.append(taxrank)
class ProteinObj():
    def __init__(self, protein_id, length, proteome_id):
        self.protein_id = protein_id
        self.proteome_id = proteome_id
        self.species_id = species_id
        self.species_name = species_name_dict.get(speciesID, None)
        self.length = length

        # OG
        self.cluster_id = None
        # gff
        self.contig_id = None
        # blast
        self.taxonomy = None # dict : key=rank, val=taxid ; translateOnDemand
        self.AI = None
        self.HI = None
        # interproscan
        self.domain_by_source = {}

        self.domain_list = None
        self.domain_set = None
        self.domain_diversity = None
        self.domain_count = None
        self.domain_counter = None

        def add_domain(self, domainObj):
            self.domain_by_source[domainObj.source]
            self.domain_list.append(domainObj.id)
            self.domain_count = len(self.domain_list)
            self.domain_counter = domain_counter.get(domainObj.id, 0) + 1

class ProteinCollection():

class DomainRecordObj(object):
    """docstring for DomainObj"""
    def __init__(self, domain_id, protein_id, domain_source, domain_evalue, domain_desc, domain_start, domain_stop):
        self.domain_id = domain_id
        self.protein_id = protein_id
        self.source = domain_source
        self.evalue = domain_evalue
        self.desc = domain_desc
        self.start = domain_start
        self.stop = domain_stop

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

    PROTEOME_IDS = {} # OrthofinderID to ProteomeID
    SEQUENCE_IDS = {} # ProteinID to ProteomeID

    try:
        species_classification_f = args['--species_classification_file']
        nodesdb_f = args['--nodesdb']
        species_ids_f = args['--species_file']
        sequence_ids_f = args['--sequence_file']
        domain_f = args['--domain_file']
    except docopt.DocoptExit:
        print __doc__.strip()

    aloCollection = AloCollection()
    aloCollection.parse_classification(species_classification_f)

    PROTEOME_IDS = parse_species_ids(species_ids_f)




    domainCollection = DomainCollection()
    print "[STATUS] - Parsing domains from %s" % (domain_f)
    domainCollection.parse_domains(domain_f)
    domainCollection.count_combinations()
