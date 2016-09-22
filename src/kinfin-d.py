#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: kinfin-d.py      [-s <FILE>] [-g <FILE>] [-c <FILE>]
                        [-d <FILE>] [-b <FILE>] [--fasta_dir <DIR>]
                        [--sequence_file <FILE>]
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
            --sequence_file <FILE>              SequenceIDs.txt used in OrthoFinder
            -g, --groups <FILE>                 OrthologousGroups.txt produced by OrthoFinder
            -c, --classification_file <FILE>    SpeciesClassification file
            -d, --interproscan_file <FILE>      InterProScan TSV file
            -b, --blast_file <FILE>             BLAST results of similarity search against NR/Uniref90
            --fasta_dir <DIR>                   Directory of FASTA files
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

class ClusterObj():
    def __init__(self, cluster_id, protein_id):
        self.cluster_id = cluster_id
        self.proteins = proteins
        try:
            self.proteomeIDs = [x.split(DELIMITER)[0] for x in proteinIDs]
        except:
            sys.exit('[ERROR] - Bad delimiter "%s"' % DELIMITER)

        self.singleton = False if len(self.proteins) > 1 else True

        self.proteinIDs_by_proteomeID = self.generate_proteinIDs_by_proteomeID()
        self.proteomeIDs_unique = set(self.proteomeIDs)
        self.proteinID_count = len(proteinIDs)
        self.proteinID_count_by_proteomeID = Counter(self.proteomeIDs)

        self.levelIDs_by_rank = {}
        self.coverage_by_levelID_by_rankID = {} # at least
        self.cluster_type_by_rankID = {}

        self.protein_length = []
        self.protein_length_median = 0
        self.filter_pass = False
        # interproscan results
        self.domain_composition_count = {}

    def generate_proteinIDs_by_proteomeID(self):
        proteinIDs_by_proteomeID = {}
        for proteinID in self.proteinIDs:
            proteomeID = proteinID.split(DELIMITER)[0]
            if not proteomeID in proteinIDs_by_proteomeID:
                proteinIDs_by_proteomeID[proteomeID] = set()
            proteinIDs_by_proteomeID[proteomeID].add(proteinID)
        return proteinIDs_by_proteomeID

class AttributeLevelObj():
    '''
    ALOs are build from config file
        - clusters can be ALO-specific or not, if they are not, they will belong to several ALOs
        - singleton clusters are always ALO-specific
        - does not save clusters but keys of clusters
    '''
    def __init__(self, attribute, level, proteomes):
        self.attribute = attribute
        self.level = level
        self.proteomes = proteomes
        self.proteome_count = len(proteomes)

        self.proteins = [] # proteins of proteomes within this ALO
        self.protein_count = 0
        self.clusters = [] # clusters containing a protein of at least one member of this ALO
        self.cluster_count = 0

        self.clusters_by_type = { 'specific' : [], 'shared' : []}

        self.clusters = {'true_1to1' : [], 'fuzzy_1to1' : []}

        self.protein_by_type = {'singleton' : [], 'monoton' : [], 'multiton' : [], 'true_1to1': {'monoton' : [], 'multiton' : []}, 'fuzzy_1to1': {'monoton' : [], 'multiton' : []}}
        self.protein_count_by_type = {'singleton' : 0, 'monoton' : 0, 'multiton' : 0, 'true_1to1' : {'monoton' : 0, 'multiton' : 0}, 'fuzzy_1to1': {'monoton' : 0, 'multiton' : 0}}
        self.cluster_by_type = {'singleton' : [], 'monoton' : [], 'multiton' : [], 'true_1to1': {'monoton' : [], 'multiton' : []}, 'fuzzy_1to1': {'monoton' : [], 'multiton' : []}}
        self.cluster_count_by_type = {'singleton' : 0, 'monoton' : 0, 'multiton' : 0, 'true_1to1' : {'monoton' : 0, 'multiton' : 0}, 'fuzzy_1to1': {'monoton' : 0, 'multiton' : 0}}
        self.protein_span = []

        self.domainIDs = []
        self.domainID_count = 0
        self.coverage_in_clusters = []

        self.rarefaction_data = {} # repetition : number of clusters

class AloCollection():
    def __init__(self):
        self.proteomes = set() # set of proteomes for orthofinder
        self.proteome_count = 0
        self.attributes = [] # list of attribute
        self.levels_by_attributes_by_proteome = {}

        self.levels_by_attribute = {}
        self.level_count_by_attribute = {}
        self.proteomes_by_level_by_attribute = {}
        self.proteome_count_by_level_by_attribute = {}

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

        # set of levels by attribute
        self.levels_by_attribute = {attribute : set() for attribute in self.attributes}
        # set of proteomes by level by attribute
        self.proteomes_by_level_by_attribute = {attribute : {} for attribute in self.attribute}
        # count of proteomes by level by attribute
        self.proteome_count_by_level_by_attribute = {attribute : {} for attribute in self.attribute}

        for proteome in self.levels_by_attributes_by_proteome:
            for attribute in self.levels_by_attributes_by_proteome[proteome]:
                for level in self.levels_by_attributes_by_proteome[proteome][attribute]:
                    self.levels_by_attribute[attribute].add(level)
                    if not level in self.proteomes_by_level_by_attribute[attribute]:
                        self.proteomes_by_level_by_attribute[attribute][level] = set()
                        self.proteome_count_by_level_by_attribute[attribute][level] = 0
                    self.proteomes_by_level_by_attribute[attribute][level].add(proteome)
                    self.proteome_count_by_level_by_attribute[attribute][level] += 1
                    self.proteome_count_by_level_by_attribute

        self.proteome_count = len(self.proteomes)
        self.level_count_by_attribute = {attribute : len(levels) for attribute, level in self.levels_by_attribute.items()}

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
    pass

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

def add_domain_combos(domain_list):
    w = 1
    l = len(domain_list) + 1
    while w < l + 1:
        for i in range(l-w):
            domain_combo = tuple(domain_list[i:i+w])
            if not domain_combo in domain_combos:
                domain_combos[domain_combo] = set()


            print "".join(list(domain_combo))
            domain_combos[domain_combo].add(set()
        w += 1
    print domain_list
    for d, c in domain_combos.items():
        print c, d



if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    PROTEOME_IDS = {} # OrthofinderID to ProteomeID
    SEQUENCE_IDS = {} # ProteinID to ProteomeID


    species_classification_f = args['--classification_file']
    nodesdb_f = args['--nodesdb']
    species_ids_f = args['--species_file'] # only needed if one wants to parse the FASTAs
    sequence_ids_f = args['--sequence_file'] # only needed if there is no prefix
    interproscan_f = args['--interproscan_file']

    #except docopt.DocoptExit:
    #    print __doc__.strip()

    generate_domain_combos(['a','a','b','c','d'])
    generate_domain_combos(['a','a'])
    generate_domain_combos(['a','a','a','a','a','a','a'])
    generate_domain_combos(['0'])

    #aloCollection = AloCollection()
    #aloCollection.parse_classification(species_classification_f)

    #PROTEOME_IDS = parse_species_ids(species_ids_f)

    #domainCollection = DomainCollection()
    #print "[STATUS] - Parsing domains from %s" % (interproscan_f)
    #domainCollection.parse_domains(interproscan_f)
    #domainCollection.count_combinations()
