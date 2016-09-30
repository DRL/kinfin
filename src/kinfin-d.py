#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: kinfin-d.py      [-s <FILE>] -g <FILE> -c <FILE>
                        [-d <FILE>] [-b <FILE>] [--fasta_dir <DIR>]
                        [-t <FILE>]
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
            -t, --tree_file <FILE>              Tree file (on which ALOs are defined)
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
import shutil

from docopt import docopt, DocoptExit
from collections import defaultdict
from collections import Counter
from itertools import chain, combinations, takewhile

from math import sqrt
#from ete2 import Tree
########################################################################
# General constructs
########################################################################

tree_dict = lambda : defaultdict(tree_dict)

########################################################################
# General functions
########################################################################

def get_lineage(taxid, NODESDB): # works
    lineage = {taxrank : 'undef' for taxrank in TAXRANKS}
    parent = ''
    node = taxid
    while not parent == "1":
        taxrank = NODESDB[node]['rank']
        name = NODESDB[node]['name']
        parent = NODESDB[node]['parent']
        if taxrank in TAXRANKS:
            lineage[taxrank] = name
        node = parent
    return lineage

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def mean(lst):
    return float(sum(lst)) / max(len(lst), 1)

def mode(lst):
    most_frequent = Counter(lst).most_common()
    modes = list(takewhile(lambda x_f: x_f[1] == most_frequent[0][1], most_frequent))
    return [x[0] if x[0] > 0 else 1 for x in modes]

def sd(lst, population=True):
    n = len(lst)
    differences = [x - mean(lst) for x in lst]
    sq_differences = [d ** 2 for d in differences]
    ssd = sum(sq_differences)
    if population is True:
        variance = ssd / n
    else:
        variance = ssd / (n - 1)
    sd = sqrt(variance)
    return sd

def is_mad_outlier(lst1, lst2):
    threshold=3.5
    total_list = lst1 + lst2
    diffs = [(x - median(total_list))**2 for x in total_list]
    diffsqrts = [sqrt(diff) for diff in diffs]
    print diffsqrts
    mad = median(diffsqrts)
    print mad
    if mad == 0:
        mad = 1
    modified_z_scores = [0.6745 * diffsqrt / mad for diffsqrt in diffsqrts]
    print modified_z_scores
    result = [True if modified_z_score >= threshold else False for modified_z_score in modified_z_scores]
    return result

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

class ClusterCollection():
    def __init__(self, cluster_id, protein_id):
        self.cluster_order = []
        self.clusters = {}

        self.cluster_count = 0
        self.singleton_count = 0

        self.protein_metrics = {'mean' : 0.0, 'mode' : 0.0, 'median' : 0.0, 'sd' : 0.0} # only non-singletons
        self.proteome_metrics = {'mean' : 0.0, 'mode' : 0.0, 'median' : 0.0, 'sd' : 0.0} # only non-singletons

class ClusterObj():
    def __init__(self, cluster_id, proteins):
        self.cluster_id = cluster_id
        self.proteins = proteins
        self.protein_count = len(proteins)
        self.singleton = False if self.protein_count > 1 else True
        try:
            self.proteomes_list = [x.split(DELIMITER)[0] for x in proteins]
        except:
            sys.exit('[ERROR] - Bad delimiter "%s"' % DELIMITER)
        self.proteomes = frozenset(self.proteomes_list)
        self.proteins_by_proteome = self.generate_proteins_by_proteome()
        self.protein_count_by_proteome = Counter(self.proteomes_list)

        self.levelIDs_by_rank = {}
        self.coverage_by_levelID_by_rankID = {} # at least
        self.cluster_type_by_rankID = {}

        self.protein_length = []
        self.protein_length_median = 0
        self.filter_pass = False
        # interproscan results
        self.domain_composition_count = {}

    def generate_proteins_by_proteome(self):
        proteins_by_proteome = {}
        for protein in self.proteins:
            proteome = protein.split(DELIMITER)[0]
            if not proteome in proteins_by_proteome:
                proteins_by_proteome[proteome] = set()
            proteins_by_proteome[proteome].add(protein)
        return proteins_by_proteome

class AttributeLevelObj():
    '''
    ALOs are build from config file
        - clusters can be ALO-specific or not, if they are not, they will belong to several ALOs
        - singleton clusters are always ALO-specific
        - does not save clusters but keys of clusters
    Definitions:
        'shared' : shared between one ALO and others
        'singleton' : cardinality of 1 ('specific', but separate)
        'specific' : only present within one ALO
        'missing'   : not present within this ALO (can be shared by others)

    '''
    def __init__(self, attribute, level, proteomes):
        self.attribute = attribute # string
        self.level = level # string
        self.proteomes = frozenset(proteomes) # frozenset(), used for checking whether cluster and ALO intersect
        self.proteome_count = len(proteomes) # int

        self.cluster_by_cluster_type = {'singleton' : [], 'specific' : [], 'shared' : [], 'missing' : []}  # sums up to cluster_count
        self.protein_count_by_cluster_type = {'singleton' : 0, 'specific' : 0, 'shared' : 0} # output as percentages of total protein count
        self.clusters_by_cluster_cardinality = {'shared' : {'true' : [], 'fuzzy' : []}, 'specific' : {'true' : [], 'fuzzy' : []}}

        self.cluster_mean_enrichment_by_cluster_id = {}
        self.cluster_mode_enrichment_by_cluster_id = {}
        self.cluster_zscore_by_cluster_id = {}

        self.protein_span = {'singleton' : 0, 'specific' : 0, 'shared' : 0}

        self.domainIDs = []
        self.domainID_count = 0

        self.coverage_in_clusters = []

        self.rarefaction_data = {} # repetition : number of clusters

    def add_clusterObj(self, clusterObj, cluster_type, cluster_cardinality, mean_enrichment,
        mode_enrichment, zscore, ALO_proteins):

        self.cluster_by_cluster_type[cluster_type].append(clusterObj.cluster_id)
        if not cluster_type == 'missing':
            self.protein_count_by_cluster_type[cluster_type] += len(ALO_proteins)
            # Protein span

        if cluster_cardinality:
            self.clusters_by_cluster_cardinality[cluster_type][cluster_cardinality].append(clusterObj.cluster_id)
        if cluster_type == 'shared':
            self.cluster_mean_enrichment_by_cluster_id[clusterObj.cluster_id] = mean_enrichment
            self.cluster_mode_enrichment_by_cluster_id[clusterObj.cluster_id] = mode_enrichment
            self.cluster_zscore_by_cluster_id[clusterObj.cluster_id] = zscore

    def get_protein_count(self):
        pass

    def get_cluster_count(self):
        pass

class AloCollection():
    def __init__(self, proteomes, proteome_id_by_species_id, attributes, level_by_attribute_by_proteome):
        self.proteomes = frozenset(proteomes) # set of proteomes for orthofinder
        self.proteome_count = len(proteomes)
        self.proteome_id_by_species_id = proteome_id_by_species_id
        self.level_by_attribute_by_proteome = level_by_attribute_by_proteome
        self.attributes = attributes # list of attributes


        self.levels_by_attribute = self.compute_levels_by_attribute()
        self.level_count_by_attribute = self.compute_level_count_by_attribute()
        self.proteomes_by_level_by_attribute = self.compute_proteomes_by_level_by_attribute()
        self.proteome_count_by_level_by_attribute = self.compute_proteome_count_by_level_by_attribute()

        self.ALO_by_level_by_attribute = self.create_ALOs()

    ###############################
    ### Setup                   ###
    ###############################

    def compute_levels_by_attribute(self):
        levels_by_attribute = {attribute : set() for attribute in self.attributes}
        for proteome in self.level_by_attribute_by_proteome:
            for attribute in self.level_by_attribute_by_proteome[proteome]:
                level = self.level_by_attribute_by_proteome[proteome][attribute]
                levels_by_attribute[attribute].add(level)
        return levels_by_attribute

    def compute_level_count_by_attribute(self):
        return {attribute : len(levels) for attribute, levels in self.levels_by_attribute.items()}

    def compute_proteomes_by_level_by_attribute(self):
        proteomes_by_level_by_attribute = {attribute : {} for attribute in self.attributes}
        for proteome in self.level_by_attribute_by_proteome:
            for attribute in self.level_by_attribute_by_proteome[proteome]:
                level = self.level_by_attribute_by_proteome[proteome][attribute]
                if not level in proteomes_by_level_by_attribute[attribute]:
                    proteomes_by_level_by_attribute[attribute][level] = set()
                proteomes_by_level_by_attribute[attribute][level].add(proteome)
        return proteomes_by_level_by_attribute

    def compute_proteome_count_by_level_by_attribute(self):
        proteome_count_by_level_by_attribute = {}
        for attribute in self.attributes:
            proteome_count_by_level_by_attribute[attribute] = {}
            for level in self.levels_by_attribute[attribute]:
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

    ###############################
    ### Functionality           ###
    ###############################

    def analyse_clusters(self, clusterCol):
        cluster_count = clusterCol.cluster_count
        print "\t Clusters found = %s" % (cluster_count)
        parse_steps = cluster_count/100
        print "[STATUS] - Analysing clusters ..."
        for idx, clusterObj in enumerate(clusterCol.clusterObjs):
            self.add_clusterObj_to_ALOs(clusterObj)
            progress(idx+1, parse_steps, cluster_count)

    def add_clusterObj_to_ALOs(self, clusterObj):
        '''This function selects the ALOs to which the cluster has to be added'''
        if clusterObj.singleton == True:
            for attribute in self.level_by_attribute_by_proteome[clusterObj.proteomes]:
                for level in self.level_by_attribute_by_proteome[clusterObj.proteomes][attribute]:
                    ALO = self.ALO_by_level_by_attribute[attribute][level]
                    ALO.add_clusterObj(clusterObj, 'singleton')
        else:
            levels_seen = set()
            levels_missing = set()
            for attribute in self.attributes:
                for proteome in clusterObj.proteomes:
                    levels_seen.add(self.level_by_attribute_by_proteome[proteome][attribute])
                for level in self.levels_by_attribute[attribute]:
                    cluster_type = None
                    cluster_cardinality = None
                    ALO_mean_enrichment_in_cluster = None
                    ALO_mode_enrichment_in_cluster = None
                    ALO_zscore_in_cluster = None
                    ALO_proteins = None
                    ALO = self.ALO_by_level_by_attribute[attribute][level]
                    if not level in levels_seen: # levels : missing
                        cluster_type = 'missing'
                    else: # levels : seen, either 'shared' or 'specific'
                        proteomes_in_ALO_in_cluster = clusterObj.proteomes.intersection(ALO.proteomes)
                        proteomes_in_ALO_in_cluster_counts = [clusterObj.protein_count_by_proteome.get(proteome, 0) for proteome in ALO.proteomes]
                        ALO_proteins = list(chain.from_iterable([proteins for proteome, proteins in clusterObj.proteins_by_proteome.items() if proteome in ALO.proteomes]))
                        #print attribute, level, ALO.proteomes
                        #print ALO_proteins

                        # cluster cardinality : true, fuzzy
                        if len(proteomes_in_ALO_in_cluster) > 1:
                            if all(count == 1 for count in proteomes_in_ALO_in_cluster_counts):
                                cluster_cardinality = 'true'
                            elif len(proteomes_in_ALO_in_cluster) >= 3:
                                proteomes_in_ALO_in_cluster_at_fuzzycount_count = len([x for x in proteomes_in_ALO_in_cluster_counts if x == FUZZY_COUNT])
                                proteomes_in_ALO_in_cluster_in_fuzzyrange_count = len([x for x in proteomes_in_ALO_in_cluster_counts if x in FUZZY_RANGE])
                                fuzzy_fraction = proteomes_in_ALO_in_cluster_at_fuzzycount_count/len(ALO.proteomes)
                                if fuzzy_fraction >= FUZZY_FRACTION:
                                    if proteomes_in_ALO_in_cluster_at_fuzzycount_count + proteomes_in_ALO_in_cluster_in_fuzzyrange_count == len(ALO.proteomes):
                                        cluster_cardinality = 'fuzzy'
                            else:
                                pass

                        if len(levels_seen) == 1: # specific
                            cluster_type = 'specific'
                        else: # shared
                            cluster_type = 'shared'
                            mean_proteomes_in_ALO_in_cluster = mean(proteomes_in_ALO_in_cluster_counts)
                            mode_proteomes_in_ALO_in_cluster = mode(proteomes_in_ALO_in_cluster_counts)
                            proteomes_not_in_ALO_in_cluster_count = [clusterObj.protein_count_by_proteome.get(proteome, 0) for proteome in self.proteomes if not proteome in ALO.proteomes]
                            mean_proteomes_not_in_ALO_in_cluster = mean(proteomes_not_in_ALO_in_cluster_count)
                            mode_proteomes_not_in_ALO_in_cluster = mode(proteomes_not_in_ALO_in_cluster_count)
                            proteomes_in_cluster_count = proteomes_in_ALO_in_cluster_counts + proteomes_not_in_ALO_in_cluster_count
                            mean_proteomes_in_cluster_count = mean(proteomes_in_cluster_count)
                            # zscore
                            ALO_zscore_in_cluster = (mean_proteomes_in_ALO_in_cluster - mean_proteomes_in_cluster_count)/sd(proteomes_in_cluster_count)
                            # mean enrichment
                            ALO_mean_enrichment_in_cluster = mean_proteomes_in_ALO_in_cluster/mean_proteomes_not_in_ALO_in_cluster if mean_proteomes_not_in_ALO_in_cluster > 0 else 1

                            # mode enrichment (if multimodal, mean of mode_enrichments)
                            # ALO_mode_enrichment_in_cluster = mean([mode_in/mode_out for mode_in in mode_proteomes_in_ALO_in_cluster for mode_out in mode_proteomes_not_in_ALO_in_cluster])

                            #if ALO_zscore_in_cluster < -0.5:
                            #    print attribute, level, ALO.proteomes
                            #    print "mean_enrich", ALO_mean_enrichment_in_cluster
                            #    print "zscore", ALO_zscore_in_cluster
                            #    print "in ALO", proteomes_in_ALO_in_cluster_counts
                            #    print "not in ALO", proteomes_not_in_ALO_in_cluster_count

                    ALO.add_clusterObj(clusterObj, cluster_type, cluster_cardinality, ALO_mean_enrichment_in_cluster,
                            ALO_mode_enrichment_in_cluster, ALO_zscore_in_cluster, ALO_proteins)

class DataFactory():
    def __init__(self):
        # Input
        self.classification_file = None
        self.species_ids_file = None
        self.sequence_ids_file = None
        self.nodesdb_file = None
        self.clusters_file = None
        self.interproscan_file = None
        self.protein_dir = None

        self.cluster_count = 0
        self.protein_count = 0
        # Output
        self.dirs = None
        self.alo_count = None
        self.alo_file = None

        self.clusterObj_file = None

    def readFastaLen(self, infile):
        with open(infile) as fh:
            header, seqs = '', []
            for l in fh:
                if l[0] == '>':
                    if (header):
                        yield header, ''.join(seqs)
                    header, seqs = l[1:-1].split()[0], [] # Header is split at first whitespace
                else:
                    seqs.append(l[:-1])
            #yield header, ''.join(seqs)
            yield header, len(seqs)

    def build_ProteinCollection(protein_dir, species_ids_f, sequence_ids_f):
        proteinObjs = []
        for line in read_file(sequence_ids_f):
            sequence_id, protein_id = line.split(": ")
            species_id = sequence_id.split("_")[0]
            proteome_id = aloCollection.proteome_id_by_species_id[species_id]
            proteinObj = ProteinObj(protein_id, proteome_id, species_id, sequence_id)
            proteinObjs.append(proteinObj)
        self.sequence_ids_file = sequence_ids_f
        proteinCollection = ProteinCollection(proteinObjs)
        print "\t Proteins found = %s" % (proteinCollection.protein_count)
        if protein_dir:
            print "[STATUS] - Parsing FASTAs ..."
            fasta_file_by_species_id = self.parse_species_ids(species_ids_f)
            fasta_len_by_protein_id = self.parse_protein_dir(protein_dir, fasta_file_by_species_id)
            print "[STATUS] - Adding FASTAs to ProteinCollection ..."
            parse_steps = proteinCollection.protein_count/100
            for idx, proteinObj in enumerate(proteinObjs.proteinObj):
                proteinObj.add_length(fasta_len_by_sequence_id[proteinObj.sequence_id])
                progress(idx+1, parse_steps, proteinCollection.protein_count)
        else:
            print "[STATUS] - No Fasta-Dir given, no AA-span information will be reported ..."
        return proteinCollection

    def parse_species_ids(self, species_ids_f):
        fasta_by_ortho_id = {}
        for line in read_file(species_ids_f):
            idx, fasta = line.split(": ")
            fasta_by_ortho_id[idx] = fasta
        self.species_ids_file = species_ids_f
        return fasta_by_ortho_id

    def parse_classification(self, species_classification_f):
        '''
        Input:
            species_classification_f : user defined CSV config file used for creating ALOs (Attribute-Level-Objects)
                - header: starts with '#'
                - first column is SpeciesID or

        Output:
            attributes : list of attributes (headers of species_classification_f)
            proteomes : set of proteomes
            level_by_attribute_by_proteome : dict of proteome => attribute => level
                - 1-to-1 relationship between proteome -> attribute -> level (each proteome has only one level for each attribute)
        '''
        self.check_file(species_classification_f)
        print "[STATUS] - Parsing SpeciesClassification file %s" % (species_classification_f)
        attributes = []
        level_by_attribute_by_proteome = {}
        proteomes = set()
        proteome_id_by_species_id = {}
        for line in read_file(species_classification_f):
            if line.startswith("#"):
                attributes = [x.strip() for x in line.lstrip("#").split(",")]
                if not 'IDX' == attributes[0] or not 'PROTEOME' == attributes[1]:
                    sys.exit("[ERROR] - First/second element have to be IDX/PROTEOME.\n\t%s" % (attributes))
            elif line.strip():
                temp = line.split(",")
                if not len(temp) == len(attributes):
                    sys.exit("[ERROR] - number of columns in line differs from header\n\t%s\n\t%s" % (attributes, temp))
                if temp[1] in proteomes:
                    sys.exit("[ERROR] - 'proteome' should be unique. %s was encountered multiple times" % (temp[0]))
                proteome = temp[1]
                proteomes.add(proteome)
                proteome_id_by_species_id[temp[0]] = temp[1]
                level_by_attribute_by_proteome[proteome] = {x : '' for x in attributes}
                for idx, level in enumerate(temp):
                    attribute = attributes[idx]
                    level_by_attribute_by_proteome[proteome][attribute] = level
            else:
                pass
        self.classification_file = species_classification_f
        return proteomes, proteome_id_by_species_id, attributes, level_by_attribute_by_proteome

    def update_classification(self, nodesdb_f, attributes, level_by_attribute_by_proteome):
        self.check_file(nodesdb_f)
        print "[+] - Parsing nodesDB %s" % (nodesdb_f)
        NODESDB = parse_nodesdb(nodesdb_f)
        for proteome in level_by_attribute_by_proteome:
            taxid = level_by_attribute_by_proteome[proteome]['taxid']
            lineage = get_lineage(taxid, NODESDB)
            # add lineage attribute/levels
            for taxrank in TAXRANKS:
                level_by_attribute_by_proteome[proteome][taxrank] = lineage[taxrank]
            # remove taxid-levels
            del level_by_attribute_by_proteome[proteome]['taxid']
        # remove taxid-attribute
        attributes.remove('taxid')
        # add taxranks to rank
        for taxrank in TAXRANKS:
            attributes.append(taxrank)
        self.nodesdb_file = nodesdb_f
        return attributes, level_by_attribute_by_proteome

    def build_AloCollection(self, species_classification_f, nodesdb_f):
        proteomes, proteome_id_by_species_id, attributes, level_by_attribute_by_proteome = self.parse_classification(species_classification_f)
        # Add taxonomy if needed
        if 'taxid' in set(attributes):
            print "[+] - Attribute 'taxid' found, inferring taxonomic ranks from nodesDB..."
            attributes, level_by_attribute_by_proteome = self.update_classification(nodesdb_f, attributes, level_by_attribute_by_proteome)
        return AloCollection(proteomes, proteome_id_by_species_id, attributes, level_by_attribute_by_proteome)

    def setup_dirs(self, out_prefix):
        self.dirs = {}
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
        for attribute in aloCollection.attributes:
            attribute_path = join(result_path, attribute)
            self.dirs[attribute] = attribute_path
            if not exists(attribute_path):
                print "\t%s" % (attribute_path)
                mkdir(attribute_path)

    def build_DomainCollection(self, interproscan_f):
        self.check_file(interproscan_f)
        domains_by_proteinID = {}
        print "[STATUS] - Parsing domains from %s" % (domain_f)
        with open(domain_f) as fh:
            for line in fh:
                temp = line.rstrip("\n").split()
                proteinID = temp[0]
                domain_type = temp[3] # CDD, PIRSF, Pfam, Phobius, ProSiteProfiles, SMART, SUPERFAMILY, SignalP_EUK, TIGRFAM, TMHMM
                domain_id = temp[4]
                evalue = temp[-3]
                stop = temp[-4]
                start = temp[-5]
                if len(" ".join(temp[5:-5])):
                    desc = "\"%s\"" % " ".join(temp[5:-5])
                else:
                    desc = None
                temp = line.rstrip("\n").split()
                if not proteinID in domains_by_proteinID:
                    domains_by_proteinID[proteinID] = set()
                domains_by_proteinID[proteinID].add(domain)
        for proteinID, domains in domains_by_proteinID.items():
            self.proteinObjs_by_proteinID[proteinID].domains = domains

    def check_file(self, infile):
        if not isfile(infile):
            sys.exit("[ERROR] - %s does not exist." % (infile))

    def build_ClusterCollection(self, groups_f):
        self.check_file(groups_f)
        print "[STATUS] - Parsing %s ... this may take a while" % (groups_f)
        clusterObjs = []
        with open(groups_f) as fh:
            for idx, line in enumerate(fh):
                temp = line.rstrip("\n").split(" ")
                clusterID, protein_string = temp[0], temp[1:]
                clusterObj = ClusterObj(clusterID, protein_string)
                clusterObjs.append(clusterObj)
                #except:
                #    sys.exit("[ERROR] - Line does not seem to contain clustered proteins\n%s") % line
        return ClusterCollection(clusterObjs)

    def parse_protein_dir(self, protein_dir, fasta_file_by_species_id):
        fasta_len_by_protein_id = {}
        for species_id, fasta_f in fasta_file_by_species_id.items():
            fasta_path = join(protein_dir, fasta_f)
            if not isfile(proteome_file):
                sys.exit("[ERROR] - %s does not exist." % (proteome_file))
            print "[STATUS] - parsing FASTA %s" % (fasta_f)
            for header, length in self.readFastaLen(fasta_path):
                fasta_len_by_protein_id[header] = length
        return fasta_len_by_protein_id

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

class ClusterCollection():
    def __init__(self, clusterObjs):
        self.clusterObjs = clusterObjs
        self.cluster_count = len(clusterObjs)

class ProteinCollection():
    def __init__(self, proteinObjs):
        self.proteinObjs = proteinObjs
        self.protein_count = len(proteinObjs)

class ProteinObj():
    def __init__(self, protein_id, proteome_id, species_id, sequence_id):
        self.protein_id = protein_id
        self.proteome_id = proteome_id
        self.species_id = species_id
        self.sequence_id = sequence_id
        self.length = None

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

    def add_length(self, length):
        self.length = length

    def add_domain(self, domainObj):
        self.domain_by_source[domainObj.source]
        self.domain_list.append(domainObj.id)
        self.domain_count = len(self.domain_list)
        self.domain_counter = domain_counter.get(domainObj.id, 0) + 1

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

    try:
        species_ids_f = args['--species_file'] # only needed if one wants to parse the FASTAs
        sequence_ids_f = args['--sequence_file'] # only needed if there is no prefix
        groups_f = args['--groups']

        species_classification_f = args['--classification_file']
        nodesdb_f = args['--nodesdb']

        interproscan_f = args['--interproscan_file']
        fasta_dir = args['--fasta_dir']

        tree_f = args['--tree_file']

        outprefix = args['--outprefix']
        delimiter = args['--delimiter']
        FUZZY_FRACTION = float(args['-f'])
        FUZZY_COUNT = int(args['-n'])
        fuzzy_min = int(args['--min'])
        fuzzy_max = int(args['--max'])
        REPETITIONS = int(args['--repetitions']) + 1
        PLOT_FORMAT = args['--plotfmt']
        FONTSIZE = int(args['--fontsize'])
        FIGSIZE = tuple(int(x) for x in args['--plotsize'].split(","))
    except docopt.DocoptExit:
        print __doc__.strip()

    FUZZY_RANGE = set([x for x in range(fuzzy_min, fuzzy_max+1) if not x == FUZZY_COUNT])
    TAXRANKS = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'superfamily', 'family', 'subfamily', 'genus', 'species']
    DELIMITER = delimiter.replace("\"", "")

    proteinCollection = None
    domainCollection = None
    aloCollection = None
    clusterCollection = None

    dataFactory = DataFactory()
    aloCollection = dataFactory.build_AloCollection(species_classification_f, nodesdb_f)

    proteinCollection = dataFactory.build_ProteinCollection(protein_dir, species_ids_f, sequence_ids_f)

    clusterCollection = dataFactory.build_ClusterCollection(groups_f)

    domainCollection = dataFactory.build_domainCollection(interproscan_f)
    proteinCollection.analyse_domains(domainCollection)

    aloCollection.analyse_clusters(clusterCollection)
    # aloCollection:
    # - for each ALO, plot rarefaction curve
    # proteinCollection: holds info for length
    # - for each ALO, plot length distribution (each on one plot)
    # - for each Attribute, plot length distribution of ALOs (together on one plot)
    # clusterCollection: holds info for clustersizes
    # - for each ALO, output cluster_counts, mean_of_mode_size
    # domainCollection: ...
    # - for each ALO
    #   - output mean_domain_diversity by cluster type

    dataFactory.setup_dirs(outprefix)
    #aloCollection = AloCollection()
    #aloCollection.parse_classification(species_classification_f)

    #PROTEOME_IDS = parse_species_ids(species_ids_f)

    #generate_domain_combos(['a','a','b','c','d'])
    #generate_domain_combos(['a','a'])
    #generate_domain_combos(['a','a','a','a','a','a','a'])
    #generate_domain_combos(['0'])

    #if len([x for x in ALO in cluster]) == speciiens in ALO:
    #    gain
    #elif len([x for x in ALO in cluster]) > speciiens in ALO:
    #    multi
    #elif len([x for x in ALO in cluster]) > speciiens in ALO:
    #else
    #domainCollection = DomainCollection()
    #print "[STATUS] - Parsing domains from %s" % (interproscan_f)
    #domainCollection.parse_domains(interproscan_f)
    #domainCollection.count_combinations()


'''
PLOT:
- do clusters with high zscores have particular domains?
- plot enrichment/zscore distribution across clusters for each ALO.
Loss of 1:1s :

    for sp in species
        - get shared clusters where it's count is 1
        -



'''
