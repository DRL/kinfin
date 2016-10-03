#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: kinfin-d.py      -p <FILE> -s <FILE> -g <FILE> -c <FILE>
                        [-d <FILE>] [-b <FILE>] [-a <DIR>]
                        [-t <FILE>]
                        [--nodesdb <FILE>] [--delimiter <STRING>]
                        [-f <FLOAT>] [-n <INT>] [--min <INT>] [--max <INT>]
                        [-l <INT>] [-r <INT>]
                        [--fontsize <INT>] [--plotsize INT,INT]
                        [-o <PREFIX>] [--plotfmt <PLOTFORMAT>]
                        [-h|--help]

    Options:
        -h --help                           show this

        Input files
            -p, --species_file <FILE>           SpeciesIDs.txt used in OrthoFinder
            -s, --sequence_file <FILE>          SequenceIDs.txt used in OrthoFinder
            -c, --attributes_file <FILE>        SpeciesClassification file
            -g, --groups_file <FILE>                 OrthologousGroups.txt produced by OrthoFinder
            -d, --interproscan_file <FILE>      InterProScan TSV file
            -b, --blast_file <FILE>             BLAST results of similarity search against NR/Uniref90
            -a, --fasta_dir <DIR>               Directory of FASTA files
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
            --plotfmt <STR>                     Plot formats [default: png]

"""

from __future__ import division
import sys
from os.path import basename, isfile, abspath, splitext, join, exists
from os import getcwd, mkdir
import shutil
import ast

from collections import Counter, Iterable, defaultdict
from itertools import chain, combinations, takewhile
from math import sqrt, log

import_errors = []
try:
    from docopt import docopt, DocoptExit
except ImportError:
    import_errors.append("[ERROR] : Module \'Docopt\' was not found. Please install \'Docopt\' using \'pip install docopt\'")
try:
    import six
except ImportError:
    import_errors.append("[ERROR] : Module \'six\' was not found. Please install \'six\' using \'pip install six\'")
try:
    from ete3 import Tree
except ImportError:
    import_errors.append("[ERROR] : Module \'ete3\'' was not found. Please install \'ete3\' using \'pip install ete3\'")
if import_errors:
    sys.exit("\n".join(import_errors))

########################################################################
# General functions
########################################################################

def check_file(infile):
        if not isfile(infile):
            sys.exit("[ERROR] - %s does not exist." % (infile))

def generate_domain_combinations(domain_list):
    w = 1
    l = len(domain_list) + 1
    while w < l + 1:
        for i in range(l-w):
            domain_combo = tuple(domain_list[i:i+w])
            if not domain_combo in domain_combos:
                domain_combos[domain_combo] = set()
            domain_combos[domain_combo].add(set())
        w += 1
    print domain_list
    return domain_combos

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

def shannon_diversity_index(count_d):
    N = sum(count_d.values())
    def p(n, N):
        if n == 0:
            return 0
        else:
            return (float(n)/N) * log(float(n)/N)
    return -sum(p(n, N) for n in count_d.values() if n is not 0)

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
        print "[PROGRESS] \t- %d%%" % (100)
    elif int(iteration) % int(steps+1) == 0:
        sys.stdout.write('\r')
        print "[PROGRESS] \t- %d%%" % (float(int(iteration)/int(max_value))*100),
        sys.stdout.flush()
    else:
        pass

def read_file(f):
    if not f or not exists(f):
        sys.exit("[ERROR] - File %s does not exist." % (f))
    with open(f) as fh:
        for line in fh:
            yield line.rstrip("\n")

########################################################################
# CLASS : DataFactory
########################################################################

class DataFactory():
    def __init__(self):
        # Input
        self.attributes_file = None
        self.species_ids_file = None
        self.sequence_ids_file = None
        self.nodesdb_file = None
        self.clusters_file = None
        self.interproscan_file = None
        self.fasta_dir = None
        self.tree_file = None

        self.cluster_count = 0
        self.protein_count = 0
        # Output
        self.dirs = None
        self.alo_count = None
        self.alo_file = None

        self.clusterObj_file = None

    ###############################
    ### build_AloCollection
    ###############################

    def build_AloCollection(self, species_attributes_f, nodesdb_f, tree_f):
        proteomes, proteome_id_by_species_id, attributes, level_by_attribute_by_proteome = self.parse_attributes(species_attributes_f)
        # Add taxonomy if needed
        if 'taxid' in set(attributes):
            print "[+] - Attribute 'taxid' found, inferring taxonomic ranks from nodesDB..."
            attributes, level_by_attribute_by_proteome = self.add_taxid_attributes(nodesdb_f, attributes, level_by_attribute_by_proteome)
        # Add ALOs from tree if provided
        if tree_f:
            tree = self.parse_tree(tree_f)
            attributes, level_by_attribute_by_proteome = self.add_attribute_levels_from_tree(tree, attributes, level_by_attribute_by_proteome)
        return AloCollection(proteomes, proteome_id_by_species_id, attributes, level_by_attribute_by_proteome)

    ###############################
    ### build_AloCollection : parse_attributes
    ###############################

    def parse_attributes(self, species_attributes_f):
        check_file(species_attributes_f)
        print "[STATUS] - Parsing SpeciesClassification file: %s ..." % (species_attributes_f)
        attributes = []
        level_by_attribute_by_proteome = {}
        proteomes = set()
        proteome_id_by_species_id = {}
        for line in read_file(species_attributes_f):
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
        self.attributes_file = species_attributes_f
        return proteomes, proteome_id_by_species_id, attributes, level_by_attribute_by_proteome

    ###############################
    ### build_AloCollection : add_taxid_attributes
    ###############################

    def add_taxid_attributes(self, nodesdb_f, attributes, level_by_attribute_by_proteome):
        check_file(nodesdb_f)
        print "[STATUS] - Parsing nodesDB %s" % (nodesdb_f)
        NODESDB = self.parse_nodesdb(nodesdb_f)
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

    ###############################
    ### build_AloCollection : add_taxid_attributes : parse_nodesdb
    ###############################

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

    ###############################
    ### build_AloCollection : parse_tree
    ###############################

    def parse_tree(self, tree_f):
        print "[STATUS] - Parsing Tree file : %s ..." % (tree_f)
        tree = Tree(tree_f)
        print tree
        #tree.render('test.pdf')
        tree = tree.write(format=9).rstrip(";")
        py_by_newick = {'(' : '[', ')' : ']', ',' : ','}
        tree_py = []
        taxon = []
        for c in tree:
            if not c in py_by_newick:
                taxon.append(c)
            else:
                if taxon:
                    tree_py.append("\'" + "".join(taxon) + "\'")
                    taxon = []
                tree_py.append(py_by_newick.get(c))
        tree = ast.literal_eval("".join(tree_py))
        print tree, isinstance(tree, list)
        return tree

    ###############################
    ### build_AloCollection : add_attribute_levels_from_tree
    ###############################

    def add_attribute_levels_from_tree(self, tree, attributes, level_by_attribute_by_proteome):
        #taxa = set([t.name for t in tree.get_leaves()])
        taxa = set([t for t in self.flatten_list(tree)])
        # check if taxa are part of PROTEOME
        proteomes = set(level_by_attribute_by_proteome.keys())
        if not taxa == proteomes:
            missing_in_tree = proteomes.difference(taxa)
            gained_in_tree = taxa.difference(proteomes)
            sys.exit("[ERROR] : Proteomes in %s differ from taxa in %s\n\tProteomes missing in tree: %s\n\tAdditional taxa in tree: %s" % (self.attributes_file, self.tree_file, missing_in_tree, gained_in_tree))
        # defining bicliques
        bicliques = [x for x in self.get_bicliques(tree, taxa, set())]
        for idx, biclique in enumerate(bicliques):
            biclique_id = "biclique_%s" % (idx)
            #print biclique_id, biclique
            attributes.append(biclique_id)
            clique, anticlique = biclique
            for proteome in clique:
                level_by_attribute_by_proteome[proteome][biclique_id] = True
            for proteome in anticlique:
                level_by_attribute_by_proteome[proteome][biclique_id] = False
        return attributes, level_by_attribute_by_proteome

    ###############################
    ### build_AloCollection : add_attribute_levels_from_tree : flatten_list
    ###############################

    def flatten_list(self, lst):
        for item in lst:
            if isinstance(item, Iterable) and not isinstance(item, basestring):
                for sub in self.flatten_list(item):
                    yield sub
            else:
                yield item

    ###############################
    ### build_AloCollection : add_attribute_levels_from_tree : get_bicliques
    ###############################

    def get_bicliques(self, tree, taxa, seen):
        for item in tree:
            biclique = []
            if isinstance(item, basestring):
                clique = set(item)
            elif isinstance(item, list):
                clique = set(self.flatten_list(item))
            else:
                print '[ERROR] : this shouldn\'t have happened...'
            anticlique = taxa.difference(clique)
            if not anticlique in seen and not clique in seen:
                seen.add(frozenset(clique))
                seen.add(frozenset(anticlique))
                biclique.append(clique)
                biclique.append(anticlique)
            if isinstance(item, list):
                for sub_clique in self.get_bicliques(item, taxa, seen):
                    if sub_clique:
                        yield sub_clique
            if biclique:
                yield biclique

    ###############################
    ### setup_dirs
    ###############################

    def setup_dirs(self, out_prefix):
        self.dirs = {}
        if (out_prefix):
            result_path = join(getcwd(), "%s.kinfin_results" % (out_prefix))
        else:
            result_path = join(getcwd(), "kinfin_results")
        self.dirs['main'] = result_path
        print "[STATUS] - Output directories in \n\t%s" % (result_path)
        if exists(result_path):
            print "[STATUS] - Directory exists. Deleting directory ..."
            shutil.rmtree(result_path)
        print "[STATUS] - Creating directories ..."
        mkdir(result_path)
        for attribute in aloCollection.attributes:
            attribute_path = join(result_path, attribute)
            self.dirs[attribute] = attribute_path
            if not exists(attribute_path):
                print "\t%s" % (attribute_path)
                mkdir(attribute_path)

    ###############################
    ### build_ProteinCollection
    ###############################

    def build_ProteinCollection(self, fasta_dir, species_ids_f, sequence_ids_f):
        proteinObjs = []
        print "[STATUS] - Parsing sequence IDs: %s ..." % sequence_ids_f
        for line in read_file(sequence_ids_f):
            sequence_id, protein_id = line.split(": ")
            species_id = sequence_id.split("_")[0]
            proteome_id = aloCollection.proteome_id_by_species_id[species_id]
            proteinObj = ProteinObj(protein_id, proteome_id, species_id, sequence_id)
            proteinObjs.append(proteinObj)
        self.sequence_ids_file = sequence_ids_f
        proteinCollection = ProteinCollection(proteinObjs)
        print "[STATUS]\t - Proteins found = %s" % (proteinCollection.protein_count)
        if fasta_dir:
            print "[STATUS] - Parsing FASTAs ..."
            fasta_file_by_species_id = self.parse_species_ids(species_ids_f)
            fasta_len_by_protein_id = self.parse_fasta_dir(fasta_dir, fasta_file_by_species_id)
            print "[STATUS] - Adding FASTAs to ProteinCollection ..."
            parse_steps = proteinCollection.protein_count/100
            for idx, proteinObj in enumerate(proteinCollection.proteinObjs):
                proteinObj.add_length(fasta_len_by_protein_id[proteinObj.protein_id])
                progress(idx+1, parse_steps, proteinCollection.protein_count)
        else:
            print "[STATUS] - No Fasta-Dir given, no AA-span information will be reported ..."
        return proteinCollection

    ###############################
    ### build_ProteinCollection : parse_species_ids
    ###############################

    def parse_species_ids(self, species_ids_f):
        fasta_by_ortho_id = {}
        for line in read_file(species_ids_f):
            idx, fasta = line.split(": ")
            fasta_by_ortho_id[idx] = fasta
        self.species_ids_file = species_ids_f
        return fasta_by_ortho_id

    ###############################
    ### build_ProteinCollection : parse_fasta_dir
    ###############################

    def parse_fasta_dir(self, fasta_dir, fasta_file_by_species_id):
        fasta_len_by_protein_id = {}
        for species_id, fasta_f in fasta_file_by_species_id.items():
            fasta_path = join(fasta_dir, fasta_f)
            if not isfile(fasta_path):
                sys.exit("[ERROR] - %s does not exist." % (fasta_path))
            print "[STATUS]\t - Parsing FASTA %s" % (fasta_path)
            for header, length in self.readFastaLen(fasta_path):
                fasta_len_by_protein_id[header] = length
        return fasta_len_by_protein_id

    ###############################
    ### build_ProteinCollection : parse_fasta_dir : readFastaLen
    ###############################

    def readFastaLen(self, infile):
        with open(infile) as fh:
            header, seqs = '', []
            for l in fh:
                if l[0] == '>':
                    if (header):
                        yield header, len(''.join(seqs))
                    header, seqs = l[1:-1].split()[0], [] # Header is split at first whitespace
                else:
                    seqs.append(l[:-1])
            #yield header, ''.join(seqs)
            yield header, len(''.join(seqs))

    ###############################
    ### build_ClusterCollection
    ###############################

    def build_ClusterCollection(self, groups_f):
        check_file(groups_f)
        print "[STATUS] - Parsing %s ... this may take a while" % (groups_f)
        clusterObjs = []
        with open(groups_f) as fh:
            for idx, line in enumerate(fh):
                temp = line.rstrip("\n").split(" ")
                clusterID, protein_string = temp[0], temp[1:]
                clusterObj = ClusterObj(clusterID, protein_string)
                clusterObjs.append(clusterObj)
        return ClusterCollection(clusterObjs)

    ###############################
    ### build_DomainCollection
    ###############################

    def build_DomainCollection(self, interproscan_f):
        check_file(interproscan_f)
        domainObjs = []
        print "[STATUS] - Parsing %s ... this may take a while" % (interproscan_f)
        with open(interproscan_f) as fh:
            for line in fh:
                temp = line.rstrip("\n").split()
                domain_source = temp[3]
                domain_id = temp[4]
                if domain_source in WHITELISTED_DOMAIN_SOURCES:
                    domain_protein_id = temp[0]
                    try:
                        domain_evalue = float(temp[-3])
                    except:
                        domain_evalue = None
                    domain_stop = int(temp[-4])
                    domain_start = int(temp[-5])
                    domain_desc = None
                    if len(" ".join(temp[5:-5])):
                        desc = "\"%s\"" % " ".join(temp[5:-5])
                    domainObj = DomainObj(domain_id, domain_source, domain_protein_id, domain_evalue, domain_start, domain_stop, domain_desc)
                    domainObjs.append(domainObj)
        return DomainCollection(domainObjs)

    ###############################
    ### write_output
    ###############################

    def write_output():
        self.write_ALO_stats()

    ###############################
    ### write_output : write_ALO_stats
    ###############################

    def write_ALO_stats():
        for attribute in aloCollection.attributes:
            ALO_stats_f = join(self.dirs[attribute], "%s.cluster_stats.txt" % (attribute))
            rank_out_string = ["\t".join([\
                "attribute", \
                "level", \
                "cluster_count", \
                "cluster_singleton_count", \
                "cluster_specific_count", \
                "cluster_shared_count", \
                "cluster_missing_count", \
                "cluster_specific_true_1to1_count", \
                "cluster_specific_fuzzy_count", \
                "cluster_shared_true_1to1_count", \
                "cluster_shared_fuzzy_count", \
                "protein_count", \
                "protein_singleton_count", \
                "protein_specific_count", \
                "protein_shared_count", \
                "protein_specific_true_1to1_count", \
                "protein_specific_fuzzy_count", \
                "protein_shared_true_1to1_count", \
                "protein_shared_fuzzy_count", \
                "protein_singleton_span", \
                "protein_specific_span", \
                "protein_shared_span", \
                "members_count", \
                "members" \
                ])]

########################################################################
# CLASS : AloCollection
########################################################################

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
    ### compute_levels_by_attribute
    ###############################

    def compute_levels_by_attribute(self):
        levels_by_attribute = {attribute : set() for attribute in self.attributes}
        for proteome in self.level_by_attribute_by_proteome:
            for attribute in self.level_by_attribute_by_proteome[proteome]:
                level = self.level_by_attribute_by_proteome[proteome][attribute]
                levels_by_attribute[attribute].add(level)
        return levels_by_attribute

    ###############################
    ### compute_level_count_by_attribute
    ###############################

    def compute_level_count_by_attribute(self):
        return {attribute : len(levels) for attribute, levels in self.levels_by_attribute.items()}

    ###############################
    ### compute_proteomes_by_level_by_attribute
    ###############################

    def compute_proteomes_by_level_by_attribute(self):
        proteomes_by_level_by_attribute = {attribute : {} for attribute in self.attributes}
        for proteome in self.level_by_attribute_by_proteome:
            for attribute in self.level_by_attribute_by_proteome[proteome]:
                level = self.level_by_attribute_by_proteome[proteome][attribute]
                if not level in proteomes_by_level_by_attribute[attribute]:
                    proteomes_by_level_by_attribute[attribute][level] = set()
                proteomes_by_level_by_attribute[attribute][level].add(proteome)
        return proteomes_by_level_by_attribute

    ###############################
    ### compute_proteome_count_by_level_by_attribute
    ###############################

    def compute_proteome_count_by_level_by_attribute(self):
        proteome_count_by_level_by_attribute = {}
        for attribute in self.attributes:
            proteome_count_by_level_by_attribute[attribute] = {}
            for level in self.levels_by_attribute[attribute]:
                proteome_count_by_level_by_attribute[attribute][level] = len(self.proteomes_by_level_by_attribute[attribute][level])
        return proteome_count_by_level_by_attribute

    ###############################
    ### create_ALOs
    ###############################

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
    ### analyse_clusters
    ###############################

    def analyse_clusters(self, clusterCol):
        cluster_count = clusterCol.cluster_count
        print "[STATUS]\t - Clusters found = %s" % (cluster_count)
        parse_steps = cluster_count/100
        print "[STATUS] - Analysing clusters ..."
        for idx, clusterObj in enumerate(clusterCol.clusterObjs):
            self.add_clusterObj_to_ALOs_and_update(clusterObj)
            progress(idx+1, parse_steps, cluster_count)

    ###############################
    ### analyse_clusters : add_clusterObj_to_ALOs_and_update
    ###############################

    def add_clusterObj_to_ALOs_and_update(self, clusterObj):
        '''This function selects the ALOs to which the cluster has to be added'''
        protein_counts_by_attribute = {}
        protein_length_mean_by_attribute = {}
        if clusterObj.singleton == True:
            for proteome_id in clusterObj.proteomes:
                for attribute in self.level_by_attribute_by_proteome[proteome_id]:
                    level = self.level_by_attribute_by_proteome[proteome_id][attribute]
                    ALO = self.ALO_by_level_by_attribute[attribute][level]
                    ALO_proteins = list(chain.from_iterable([proteins for proteome, proteins in clusterObj.protein_ids_by_proteome.items() if proteome in ALO.proteomes]))
                    ALO_protein_spans = [proteinCollection.proteinObjs_by_protein_id[protein_id].length for protein_id in ALO_proteins]
                    if not attribute in protein_counts_by_attribute:
                        protein_counts_by_attribute[attribute] = []
                        protein_length_mean_by_attribute[attribute] = []
                    protein_counts_by_attribute[attribute].append(len(ALO_proteins))
                    protein_length_mean_by_attribute[attribute].append(mean(ALO_protein_spans))
                    ALO.add_clusterObj(clusterObj, 'singleton', None, None, None, ALO_proteins, ALO_protein_spans)
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
                    ALO_protein_spans = None
                    ALO = self.ALO_by_level_by_attribute[attribute][level]
                    if not level in levels_seen: # levels : missing
                        cluster_type = 'missing'
                    else: # levels : seen, either 'shared' or 'specific'
                        proteomes_in_ALO_in_cluster = clusterObj.proteomes.intersection(ALO.proteomes)
                        proteomes_in_ALO_in_cluster_counts = [clusterObj.protein_count_by_proteome.get(proteome, 0) for proteome in ALO.proteomes]
                        ALO_proteins = list(chain.from_iterable([proteins for proteome, proteins in clusterObj.protein_ids_by_proteome.items() if proteome in ALO.proteomes]))
                        ALO_protein_spans = [proteinCollection.proteinObjs_by_protein_id[protein_id].length for protein_id in ALO_proteins]
                        if not attribute in protein_counts_by_attribute:
                            protein_counts_by_attribute[attribute] = []
                            protein_length_mean_by_attribute[attribute] = []
                        protein_counts_by_attribute[attribute].append(len(ALO_proteins))
                        protein_length_mean_by_attribute[attribute].append(mean(ALO_protein_spans))
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
                            stdev = sd(proteomes_in_cluster_count)
                            if stdev == 0:
                                ALO_zscore_in_cluster = 1.0
                            else:
                                ALO_zscore_in_cluster = (mean_proteomes_in_ALO_in_cluster - mean_proteomes_in_cluster_count)/sd(proteomes_in_cluster_count)
                            # mean enrichment
                            ALO_mean_enrichment_in_cluster = mean_proteomes_in_ALO_in_cluster/mean_proteomes_not_in_ALO_in_cluster if mean_proteomes_not_in_ALO_in_cluster > 0 else 1

                    ALO.add_clusterObj(clusterObj, cluster_type, cluster_cardinality, ALO_mean_enrichment_in_cluster, ALO_zscore_in_cluster, ALO_proteins, ALO_protein_spans)
        clusterObj.update_clusterObj(protein_counts_by_attribute, protein_length_mean_by_attribute)
        return clusterObj

########################################################################
# CLASS : AttributeLevelObj
########################################################################

class AttributeLevelObj():
    '''
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

        self.domain_count_by_domain_id_by_domain_source = None
        self.coverage_in_clusters = []
        self.rarefaction_data = {} # repetition : number of clusters

    ###############################
    ### add_clusterObj
    ###############################

    def add_clusterObj(self, clusterObj, cluster_type, cluster_cardinality, mean_enrichment, zscore, ALO_proteins, ALO_protein_spans):
        self.cluster_by_cluster_type[cluster_type].append(clusterObj.cluster_id)
        if not cluster_type == 'missing':
            self.protein_count_by_cluster_type[cluster_type] += len(ALO_proteins)
            # Protein span
            self.protein_span[cluster_type] += sum(ALO_protein_spans)
        if cluster_cardinality:
            self.clusters_by_cluster_cardinality[cluster_type][cluster_cardinality].append(clusterObj.cluster_id)
        if cluster_type == 'shared':
            self.cluster_mean_enrichment_by_cluster_id[clusterObj.cluster_id] = mean_enrichment
            self.cluster_zscore_by_cluster_id[clusterObj.cluster_id] = zscore

    ###############################
    ### get_protein_count
    ###############################

    def get_protein_count(self):
        pass

    ###############################
    ### get_cluster_count
    ###############################

    def get_cluster_count(self):
        pass

########################################################################
# CLASS : ProteinCollection
########################################################################

class ProteinCollection():
    def __init__(self, proteinObjs):
        self.proteinObjs = proteinObjs
        self.proteinObjs_by_protein_id = {proteinObj.protein_id : proteinObj for proteinObj in proteinObjs}
        self.protein_count = len(proteinObjs)
        self.domain_sources = None

    ###############################
    ### analyse_domains
    ###############################

    def analyse_domains(self, domainCollection):
        self.domain_sources = domainCollection.domain_sources
        for domainObj in domainCollection.domainObjs:
            proteinObj = self.proteinObjs_by_protein_id[domainObj.domain_protein_id]
            proteinObj.add_domainObj(domainObj)
        for proteinObj in self.proteinObjs:
            proteinObj.analyse_domainObjs()
            self.count_domains(proteinObj)

    ###############################
    ### count_domains
    ###############################

    def count_domains(self, proteinObj):
        pass

########################################################################
# CLASS : ProteinObj
########################################################################

class ProteinObj():
    def __init__(self, protein_id, proteome_id, species_id, sequence_id):
        self.protein_id = protein_id
        self.proteome_id = proteome_id
        self.species_id = species_id
        self.sequence_id = sequence_id
        self.length = None

        self.cluster_id = None
        # gff
        self.contig_id = None
        # blast
        self.taxonomy = None # dict : key=rank, val=taxid ; translateOnDemand
        self.AI = None
        self.HI = None

        # interproscan
        self.domain_list = None
        self.secreted = False

        self.domain_count_by_domain_source = None
        self.domain_count_by_domain_id_by_domain_source = None
        self.domain_diversity_by_domain_source = None

    ###############################
    ### add_length
    ###############################

    def add_length(self, length):
        self.length = length

    ###############################
    ### get_domain_list
    ###############################

    def get_domain_list(self):
        return sorted(self.domain_list, key=lambda x: x.domain_start, reverse=False)

    ###############################
    ### add_domainObj
    ###############################

    def add_domainObj(self, domainObj):
        if not self.domain_list:
            self.domain_list = []
        self.domain_list.append(domainObj)
        if domainObj.domain_id == 'SignalP-noTM':
            self.secreted = True

    ###############################
    ### analyse_domainObjs
    ###############################

    def analyse_domainObjs(self):
        domains_by_domain_source = {}
        if self.domain_list:
            for domainObj in self.domain_list:
                if not domainObj.domain_source in domains_by_domain_source:
                    domains_by_domain_source[domainObj.domain_source] = []
                domains_by_domain_source[domainObj.domain_source].append(domainObj.domain_id)
            domain_count_by_domain_source = {}
            domain_count_by_domain_id_by_domain_source = {}
            domain_diversity_by_domain_source = {}
            for source in domains_by_domain_source:
                domain_count_by_domain_source[source] = len(domains_by_domain_source[source])
                domain_count_by_domain_id_by_domain_source[source] = Counter(domains_by_domain_source[source])
                domain_diversity_by_domain_source[source] = shannon_diversity_index(domain_count_by_domain_id_by_domain_source[source])
            self.domain_count_by_domain_source = domain_count_by_domain_source
            self.domain_count_by_domain_id_by_domain_source = domain_count_by_domain_id_by_domain_source
            self.domain_diversity_by_domain_source = domain_diversity_by_domain_source

########################################################################
# CLASS : ClusterCollection
########################################################################

class ClusterCollection():
    def __init__(self, clusterObjs):
        self.clusterObjs = clusterObjs
        self.cluster_count = len(clusterObjs)
        self.protein_stats_by_attribute = None

    ###############################
    ### compute_proteins_stats_by_attribute
    ###############################

    def compute_proteins_stats_by_attribute(self):
        ''' compute overall mean/median/sd of counts and mean/sd of length across all clusterObjs'''
        pass

########################################################################
# CLASS : ClusterObj
########################################################################

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
        self.protein_ids_by_proteome = self.compute_protein_ids_by_proteome()
        self.protein_count_by_proteome = Counter(self.proteomes_list)
        self.protein_stats_by_attribute = None

        self.coverage_by_levelID_by_rankID = {} # at least
        self.domain_composition_count = {}

    ###############################
    ### compute_protein_ids_by_proteome
    ###############################

    def compute_protein_ids_by_proteome(self):
        ''' this has to be fixed to work only with species_ids, sequence_ids '''
        proteins_by_proteome = {}
        for protein in self.proteins:
            proteome = protein.split(DELIMITER)[0]
            if not proteome in proteins_by_proteome:
                proteins_by_proteome[proteome] = set()
            proteins_by_proteome[proteome].add(protein)
        return proteins_by_proteome

    ###############################
    ### update_clusterObj
    ###############################

    def update_clusterObj(self, protein_counts_by_attribute, protein_length_mean_by_attribute):
        protein_stats_by_attribute = {}
        for attribute in protein_stats_by_attribute:
            protein_stats_by_attribute[attribute] = {}
            protein_stats_by_attribute[attribute]['mean_count'] = mean(protein_stats_by_attribute[attribute])
            protein_stats_by_attribute[attribute]['median_count'] = median(protein_stats_by_attribute[attribute])
            protein_stats_by_attribute[attribute]['sd_count'] = sd(protein_stats_by_attribute[attribute])
            protein_stats_by_attribute[attribute]['mean_length'] = mean(protein_stats_by_attribute[attribute])
            protein_stats_by_attribute[attribute]['median_length'] = median(protein_stats_by_attribute[attribute])
            protein_stats_by_attribute[attribute]['sd_length'] = sd(protein_stats_by_attribute[attribute])
        self.protein_stats_by_attribute = protein_stats_by_attribute

########################################################################
# CLASS : DomainCollection
########################################################################

class DomainCollection():
    def __init__(self, domainObjs):
        self.domainObjs = domainObjs
        self.domain_count = len(domainObjs)
        self.domain_sources = set([domainObj.domain_source for domainObj in self.domainObjs])

    def compute_domain_combination_counts(self):
        pass

########################################################################
# CLASS : DomainObj
########################################################################

class DomainObj():
    def __init__(self, domain_id, domain_source, domain_protein_id, domain_evalue, domain_start, domain_stop, domain_desc):
        self.domain_id = domain_id
        self.domain_source = domain_source
        self.domain_protein_id = domain_protein_id
        self.domain_evalue = domain_evalue
        self.domain_start = domain_start
        self.domain_stop = domain_stop
        self.domain_length = domain_stop - domain_start
        self.domain_desc = domain_desc

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    try:
        species_ids_f = args['--species_file'] # only needed if one wants to parse the FASTAs
        sequence_ids_f = args['--sequence_file'] # only needed if there is no prefix
        groups_f = args['--groups_file']

        species_attributes_f = args['--attributes_file']
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
    except:
        print __doc__.strip()

    FUZZY_RANGE = set([x for x in range(fuzzy_min, fuzzy_max+1) if not x == FUZZY_COUNT])
    TAXRANKS = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'superfamily', 'family', 'subfamily', 'genus', 'species']
    WHITELISTED_DOMAIN_SOURCES = set(['Pfam', 'CDD', 'SignalP_EUK', 'SUPERFAMILY', 'SMART', 'ProSiteProfiles'])
    DELIMITER = delimiter.replace("\"", "")

    aloCollection = None
    proteinCollection = None
    domainCollection = None
    clusterCollection = None

    dataFactory = DataFactory()
    aloCollection = dataFactory.build_AloCollection(species_attributes_f, nodesdb_f, tree_f)
    dataFactory.setup_dirs(outprefix)
    proteinCollection = dataFactory.build_ProteinCollection(fasta_dir, species_ids_f, sequence_ids_f)
    if interproscan_f:
        domainCollection = dataFactory.build_DomainCollection(interproscan_f)
        proteinCollection.analyse_domains(domainCollection)
    clusterCollection = dataFactory.build_ClusterCollection(groups_f)
    aloCollection.analyse_clusters(clusterCollection)
    dataFactory.write_output()
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

'''
    #aloCollection = AloCollection()
    #aloCollection.parse_attributes(species_attributes_f)

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

PLOT:
- do clusters with high zscores have particular domains?
- plot enrichment/zscore distribution across clusters for each ALO.
Loss of 1:1s :

    for sp in species
        - get shared clusters where it's count is 1
        -



'''
