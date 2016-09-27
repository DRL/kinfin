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
from statistics import mean, mode
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

<<<<<<< HEAD
class ClusterCollection():
    def __init__(self, cluster_id, protein_id):
        self.cluster_order = []
        self.clusters = {}

        self.cluster_count = 0
        self.singleton_count = 0

        self.protein_metrics = {'mean' : 0.0, 'mode' : 0.0, 'median' : 0.0, 'sd' = 0.0} # only non-singletons
        self.proteome_metrics = {'mean' : 0.0, 'mode' : 0.0, 'median' : 0.0, 'sd' = 0.0} # only non-singletons

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
        self.attribute = attribute # string
        self.level = level # string
        self.proteomes = proteomes # set()
        self.proteome_count = len(proteomes) # int

        self.clusters_by_cluster_type = {'singleton' : [], 'gained' : [], 'shared' : [], 'lost' : []}  # sums up to cluster_count
        self.proteins_by_cluster_type = {'singleton' : [], 'gained' : [], 'shared' : [], 'lost' : []}

        self.clusters_by_cluster_cardinality = {'shared' : {'true_1_to_1' : [], 'fuzzy_n_to_n' : [], 'expanded' : [], 'contracted' : []}, 'gained' : {'true_1_to_1' : [], 'fuzzy_n_to_n' : []} # expanded/contracted only applies to shared

        self.protein_span = []

        self.domainIDs = []
        self.domainID_count = 0

        self.coverage_in_clusters = []

        self.rarefaction_data = {} # repetition : number of clusters

    def add_cluster(self, clusterObj):
        '''
            - cluster gets added to each ALO depending whether members are shared : pool of ALOs
            - this gets checked one clustertype is determined.

            if not cluster is singleton:

            else
                ALO
        '''

    def get_protein_count(self):
        pass

    def get_cluster_count(self, arg):
        if arg == 'singleton':
            pass

class AloCollection():
    def __init__(self, proteomes, attributes, levels_by_attributes_by_proteome):
        self.proteomes = proteomes # set of proteomes for orthofinder
        self.proteome_count = len(proteomes)
        self.levels_by_attributes_by_proteome = levels_by_attributes_by_proteome
        self.attributes = attributes # list of attributes
        self.levels = list(set(levels_by_attributes_by_proteome.values()))

        self.levels_by_attribute = self.compute_levels_by_attribute()
        self.level_count_by_attribute = self.compute_level_count_by_attribute()
        self.proteomes_by_level_by_attribute = self.compute_proteomes_by_level_by_attribute()
        self.proteome_count_by_level_by_attribute = self.compute_proteome_count_by_level_by_attribute()

        self.ALO_by_level_by_attribute = self.create_ALOs()

    ### Setup ###

    def compute_levels_by_attribute(self):
        levels_by_attribute = {attribute : set() for attribute in self.attributes}
        for proteome in self.levels_by_attributes_by_proteome:
            for attribute in self.levels_by_attributes_by_proteome[proteome]:
                for level in levels_by_attributes_by_proteome[proteome][attribute]:
                    levels_by_attribute[attribute].add(level)
        return levels_by_attribute

    def compute_level_count_by_attribute(self):
        return {attribute : len(levels) for attribute, levels in self.levels_by_attribute.items()}

    def compute_proteomes_by_level_by_attribute(self):
        proteomes_by_level_by_attribute = {attribute : {} for attribute in self.attributes}
        for proteome in self.levels_by_attributes_by_proteome:
            for attribute in self.levels_by_attributes_by_proteome[proteome]:
                for level in self.levels_by_attributes_by_proteome[proteome][attribute]:
                    if not level in proteomes_by_level_by_attribute[attribute]:
                        proteomes_by_level_by_attribute[attribute][level] = set()
                    proteomes_by_level_by_attribute[attribute][level].add(proteome)
        return proteomes_by_level_by_attribute

    def compute_proteome_count_by_level_by_attribute(self):
        proteome_count_by_level_by_attribute = {}
        for attribute in self.attributes:
            proteome_count_by_level_by_attribute[attribute] = {}
            for level in self.levels_by_attribute:
                proteome_count_by_level_by_attribute[attribute][level] = len(self.proteomes_by_level_by_attribute[attribute][level])
        return proteome_count_by_level_by_attribute

    def create_ALOs(self):
        ALO_by_level_by_attribute = {attribute: {} for attribute in self.attributes}
        for attribute in self.attributes:
            for level in self.levels_by_attribute[attribute]:
                proteomes = self.proteomes_by_level_by_attribute[attribute][level]
                ALO = AttributeLevelObj(attribute, level, proteomes)
                if not level in ALO_by_level_by_attribute[attribute]:
                    ALO_by_level_by_attribute[attribute][level] = {}
                ALO_by_level_by_attribute[attribute][level] = ALO
        return ALO_by_level_by_attribute

class DataFactory():
    def __init__(self):
        # Input
        self.classification_file = None
        self.nodesdb_file = None
        self.clusters_file = None
        # Output
        self.dirs = None

        self.alo_count = None
        self.alo_file = None

        self.clusterObj_count = 0
        self.clusterObj_file = None

    def parse_classification(self):
        '''
        species_classification_f : user defined CSV config file used for creating ALOs (Attribute-Level-Objects)
            - header: starts with '#', required fields: 'proteome_idx' (used for linking Orthofinder-Proteome/ProteinIDs)
        '''
        print "[STATUS] - Parsing SpeciesClassification file %s" % (species_classification_f)
        attributes = []
        levels_by_attributes_by_proteome = {}
        proteomes = set()
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
                proteomes.add(proteome)
                levels_by_attributes_by_proteome[proteome] = {x : '' for x in attributes}
                for idx, level in enumerate(temp):
                    attribute = attributes[idx]
                    levels_by_attributes_by_proteome[proteome][attribute] = level
            else:
                pass
        self.classification_file = species_classification_f
        return proteomes, attributes, levels_by_attributes_by_proteome

    def update_classification(self, attributes, levels_by_attributes_by_proteome):
        if not nodesdb_f:
            sys.exit("[ERROR] - Please specify a nodesDB file, so that taxonomic ranks can be inferred based on taxIDs\n")
        print "[+] - Parsing nodesDB %s" % (nodesdb_f)
        NODESDB = parse_nodesdb(nodesdb_f)
        for proteome in levels_by_attributes_by_proteome:
            taxid = levels_by_attributes_by_proteome[proteome]['taxid']
            lineage = get_lineage(taxid)
            # add lineage attribute/levels
            for taxrank in TAXRANKS:
                levels_by_attributes_by_proteome[proteome][taxrank] = lineage[taxrank]
            # remove taxid-levels
            del levels_by_attributes_by_proteome[proteome]['taxid']
        # remove taxid-attribute
        attributes.remove('taxid')
        # add taxranks to rank
        for taxrank in TAXRANKS:
            attributes.append(taxrank)
        self.nodesdb_file = nodesdb_f
        return attributes, levels_by_attributes_by_proteome

    def build_AloCollection(self, species_classification_f):
        proteomes, attributes, levels_by_attributes_by_proteome = self.parse_classification(species_classification_f)
        # Add taxonomy if needed
        if 'taxid' in set(attributes):
            print "[+] - Attribute 'taxid' found, inferring taxonomic ranks from nodesDB..."
            attributes, levels_by_attributes_by_proteome = self.update_classification()
        return AloCollection(proteomes, attributes, levels_by_attributes_by_proteome)

    def setup_dirs(self, out_prefix):
        self.dirs = {}
        result_path = ''
        if (out_prefix):
            result_path = join(getcwd(), "%s.kinfin_results" % (out_prefix))
        else:
            result_path = join(getcwd(), "kinfin_results")
        self.dirs['main'] = result_path
        print "[STATUS] - Output directories in \n\t%s" % (result_path)
        if exists(result_path):
            print "[STATUS] - Directory exists. Deleting directory"
            shutil.rmtree(result_path)
        print "[STATUS] - Creating directory"
        mkdir(result_path)
        for rankID in self.rankIDs:
            rank_path = join(result_path, rankID)
            self.dirs[rankID] = rank_path
            if not exists(rank_path):
                print "\t%s" % (rank_path)
                mkdir(rank_path)

    def parse_clusters(self, groups_f):
        if not isfile(groups_f) and groups_f.endswith(".txt"):
            sys.exit("[ERROR] - %s does not exist." % (groups_f))
        print "[STATUS] - Parsing %s" % (groups_f)
        clusterObjs = []
        with open(groups_f) as fh:
            for idx, line in enumerate(fh):
                try:
                    temp = line.rstrip("\n").split(" ")
                    clusterID, protein_string = temp[0], temp[1:]
                    clusterObj = ClusterObj(clusterID, protein_string)
                    clusterObjs.append(clusterObj)
                except:
                    sys.exit("[ERROR] - Line does not seem to contain clustered proteins\n%s") % line
        return ClusterCollection(clusterObjs)

    def analyse_clusters(self):
        if self.clusterObj_count:
            clusterObj_count = clusterCollection.clusterObj_count
            print "\t Clusters found = %s" % (clusterObj_count)
            parse_steps = clusterObj_count/100
            print "[STATUS] - Analysing clusters ..."
            for idx, clusterObj in enumerate(clusterCollection.clusterObjs):
                self.add_clusterObj(clusterObj)
                progress(idx+1, parse_steps, clusterObj_count)

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

class ClusterCollection():
    def __init__(self, clusterObjs):
        self.clusterObjs = clusterObjs
        self.clusterObj_count = len(clusterObjs)

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
            domain_combos[domain_combo].add(set())
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

    if len([x for x in ALO in cluster]) == speciiens in ALO:
        gain
    elif len([x for x in ALO in cluster]) > speciiens in ALO:
        multi
    elif len([x for x in ALO in cluster]) > speciiens in ALO:
    else

    dataFactory = DataFactory()
    dataFactory.setup_dirs(outprefix)
    aloCollection = dataFactory.build_AloCollection(species_classification_f)
    clusterCollection = dataFactory.parse_clusters(groups_f)
    #aloCollection = AloCollection()
    #aloCollection.parse_classification(species_classification_f)

    #PROTEOME_IDS = parse_species_ids(species_ids_f)

    #domainCollection = DomainCollection()
    #print "[STATUS] - Parsing domains from %s" % (interproscan_f)
    #domainCollection.parse_domains(interproscan_f)
    #domainCollection.count_combinations()
