#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: kb.py        -s <SPECIESIDFILE> -c <SPECIESCLASSFILE> -g <GROUPSFILE>
                    [-d FASTADIR] [-l PROTEINLENTHRESHOLD] [-r REPETITIONS]
                    [-o PREFIX] [-p PLOTFORMAT]
                    [-h|--help]

    Options:
        -h --help                           show this
        <TXT>...                            KMC2 histo txt file(s).
        --min <MIN>                         Minimum k-mer multiplicity (x) [default: 3]
        --max <MAX>                         Maximum k-mer multiplicity (x) [default: 500]
        -o, --out <OUTFILE>                 Output prefix [default: kmer_histogram]

"""

REPETITIONS = 30 + 1
PLOT_FORMAT = 'png'
species_ids_f, species_classification_f, groups_f = '','',''
try:
    species_ids_f = sys.argv[1]
    species_classification_f = sys.argv[2]
    groups_f = sys.argv[3]
from __future__ import division
import sys
from os.path import basename, isfile, abspath, splitext, join, exists
from os import getcwd, mkdir
import numpy as np
import itertools
import operator
import random
from collections import Counter
import collections

import matplotlib as mat
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
mat.use('agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import seaborn as sns
sns.set_context("talk")
#sns.set(font='serif')
sns.set(style="whitegrid")

sys.setrecursionlimit(10000) # needed for clustering
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from scipy import arange

'''
Improvements
- randomness
    - randomise proteomeID for each protein in clusters
    - plot as greys in coverage-decay plot
    (-randomise RLO membership of proteomeIDs in cluster)

'''

def readFasta(infile):
    if not isfile(infile):
         BtLog.error('0', infile)
    with open(infile) as fh:
        header, seqs = '', []
        for l in fh:
            if l[0] == '>':
                if (header):
                    yield header, ''.join(seqs)
                header, seqs = l[1:-1].split()[0], [] # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, ''.join(seqs)

############################################################################################
# DATA
############################################################################################

class DataObj():
    def __init__(self):
        # Categories
        self.rankIDs = [] # list of rankIDs
        self.rankIDs_count = 0
        self.proteomeIDs = [] # list of proteomeIDs
        self.proteomeIDs_count = 0
        self.levelIDs_by_rankID = {} # key=rankID, values=set(levelIDs)
        self.levelIDs_count_by_rankID = {} # key=rank, values=len(set(levelIDs))
        self.levelIDs_by_rankID_by_proteomeID = {} # key1=proteomeID, key2=rankID, value=levelID
        self.proteomeIDs_by_levelID_by_rankID = {} # key1=rankID, key2=levelID, value=set(proteomeIDs)
        self.count_by_levelID_by_rankID = {}
        # RankLevelObjs
        self.RLO_by_levelID_by_rankID = {} # key1=rank, key2=level, value=RankLevelObj

        # ProteinObjs
        self.proteinObjs_by_proteinID = {}

        # clusterObj
        self.clusterObjs_by_clusterID = {}
        self.clusterObjs_order = []

        # heatmap

        # misc
        self.inflation_value = ''
        self.dirs = {}
        self.fasta_parsed = False

    ############################################################################################
    # Parsing Input files
    ############################################################################################

    def parse_species_classification(self, species_classification_f):
        with open(species_classification_f) as fh:
            for l in fh:
                if l.startswith("#"):
                    self.rankIDs = [x.strip() for x in l.lstrip("#").rstrip("\n").split()]
                    self.rankIDs_count = len(self.rankIDs)
                    self.levelIDs_by_rankID = {rankID : set() for rankID in self.rankIDs}
                    self.proteomeIDs_by_levelID_by_rankID = {rankID : {} for rankID in self.rankIDs}
                else:
                    line = l.rstrip("\n").split()
                    proteomeID = line[0]
                    self.proteomeIDs.append(proteomeID)
                    self.levelIDs_by_rankID_by_proteomeID[proteomeID] = {x : '' for x in self.rankIDs}
                    for idx, levelID in enumerate(line):
                        rankID = self.rankIDs[idx]
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

    def parse_species_ids(self, species_ids_f):
        # get files for proteomes for parsing of proteins
        proteome_files = []
        with open(species_ids_f) as fh:
            for l in fh:
                number, species_fasta = l.rstrip("\n").split(": ")
                proteome_files.append(species_fasta)
        if not self.proteomeIDs_count == len(proteome_files):
            sys.exit("[ERROR] - %s fasta files found in %s. Should be %s." % (len(proteome_files), species_ids_f, self.proteomeIDs_count) )
        for proteome_file, proteome_RLO in zip(proteome_files, self.yield_RLOs(ranks=['proteome'], levels=['all'])):
                proteome_RLO.file = proteome_file

    def parse_fasta(self):
        for proteome_RLO in self.yield_RLOs(ranks=['proteome'], levels=['all']):
            if not isfile(proteome_RLO.file):
                print "[WARN] - %s does not exist." % (proteome_RLO.file)
            else:
                for header, sequence in readFasta(proteome_RLO.file):
                    proteinObj = ProteinObj(header, len(sequence), proteome_RLO.levelID)
                    self.proteinObjs_by_proteinID[proteinObj.proteinID] = proteinObj
                self.fasta_parsed = True

    def parse_clusters(self, groups_f):
        if not isfile(groups_f) and groups_f.endswith(".txt"):
            sys.exit("[ERROR] - %s does not exist." % (groups_f))
        inflation_value = basename(groups_f).lstrip("OrthologousGroups_").rstrip(".txt")
        print "[STATUS] - Parsing %s : inflation value = %s" % (groups_f, inflation_value)
        with open(groups_f) as fh:
            for line in fh:
                clusterID, protein_string = line.rstrip("\n").split(": ")
                clusterObj = ClusterObj(clusterID, protein_string.split())
                self.add_clusterObj(clusterObj)

    ############################################################################################
    # Functionality
    ############################################################################################

    def create_RLOs(self):
        for rankID in self.rankIDs:
            for levelID in self.levelIDs_by_rankID[rankID]:
                proteomeIDs = self.proteomeIDs_by_levelID_by_rankID[rankID][levelID]
                rankLevelObj = RankLevelObj(levelID, rankID, proteomeIDs)
                if not rankID in self.RLO_by_levelID_by_rankID:
                    self.RLO_by_levelID_by_rankID[rankID] = {}
                self.RLO_by_levelID_by_rankID[rankID][levelID] = rankLevelObj

    def yield_RLOs(self, **kwargs):
        rankIDs_arg = kwargs['ranks']
        #print "rankIDs_arg", rankIDs_arg
        levelIDs_arg = kwargs['levels']
        #print "levelIDs_arg", levelIDs_arg
        rankIDs = []
        if rankIDs_arg == ['all']:
            rankIDs = self.rankIDs
        else:
            for rankID in rankIDs_arg:
                if not rankID in self.rankIDs:
                    sys.exit("[ERROR] - '%s' not in ranks : %s" % (rankID, ", ".join(self.rankIDs)))
                else:
                    rankIDs.append(rankID)
        #print "rankIDs", rankIDs
        for rankID in rankIDs:
            levelIDs = []
            if levelIDs_arg == ['all']:
                levelIDs = self.levelIDs_by_rankID[rankID]
            else:
                for levelID in levelIDs_arg:
                    if not levelID in self.levelIDs_by_rankID[rankID]:
                        levelIDs_at_rank = ", ".join(list(self.levelIDs_by_rankID[rankID]))
                        sys.exit("[ERROR] - '%s' not in levels : %s at rank" % (levelID, levelIDs_at_rank, rankID))
                    else:
                        levelIDs.append(levelID)
            #print "levelIDs", levelIDs
            for levelID in levelIDs:
                yield self.RLO_by_levelID_by_rankID[rankID][levelID]

    def yield_clusterObj(self):
        for clusterID in self.clusterObjs_order:
            yield self.clusterObjs_by_clusterID[clusterID]

    def add_clusterObj(self, clusterObj):
        self.clusterObjs_order.append(clusterObj.clusterID)
        self.clusterObjs_by_clusterID[clusterObj.clusterID] = clusterObj
        # determine length parameters
        if self.fasta_parsed == True:
            for proteinID in clusterObj.proteinIDs:
                protein_length = self.proteinObjs_by_proteinID[proteinID].length
                clusterObj.protein_length.append(protein_length)
            clusterObj.protein_length_median = np.median(clusterObj.protein_length)
            if clusterObj.protein_length_median >= MEDIAN_LENGTH_THRESHOLD:
                clusterObj.filter_pass = True
        else:
            clusterObj.filter_pass = True
        # dealing with cluster
        if clusterObj.filter_pass == True:
            for rankID in self.rankIDs:
                # determine clustertype for each rank
                levelIDs_seen = set()
                for proteomeID in clusterObj.proteomeIDs:
                    levelIDs_seen.add(self.levelIDs_by_rankID_by_proteomeID[proteomeID][rankID])
                if clusterObj.proteinID_count == 1:
                    clusterObj.cluster_type_by_rankID[rankID] = 'singleton'
                else:
                    if len(levelIDs_seen) > 1:
                        clusterObj.cluster_type_by_rankID[rankID] = 'multiton'
                    else:
                        clusterObj.cluster_type_by_rankID[rankID] = 'monoton'

                for levelID in self.levelIDs_by_rankID[rankID]:
                    #if not rankID in self.protein_count_by_levelID_by_rankID:
                    #    self.protein_count_by_levelID_by_rankID[rankID] = []
                    # this has to be done by levelID not by proteomeID
#                    self.protein_count_by_levelID_by_rankID[rankID].append(clusterObj.proteinID_count_by_proteomeID)
                    if levelID in levelIDs_seen:
                        # update count in RLOs
                        for RLO in self.yield_RLOs(ranks=[rankID], levels=[levelID]):
                            # coverage
                            if not clusterObj.cluster_type_by_rankID[rankID] == 'singleton':
                                proteomeIDs_in_cluster = clusterObj.proteomeIDs_unique
                                proteomeIDs_in_RLO = RLO.proteomeIDs
                                RLO_proteomeIDs_in_cluster = proteomeIDs_in_cluster.intersection(proteomeIDs_in_RLO)
                                coverage = float(len(RLO_proteomeIDs_in_cluster)/len(proteomeIDs_in_RLO))
                                RLO.coverage_in_clusters.append(coverage)
                            # clusters
                            RLO.clusterIDs.append(clusterObj.clusterID)
                            RLO.clusterID_count += 1
                            cluster_type = clusterObj.cluster_type_by_rankID[rankID]
                            RLO.clusterID_by_type[cluster_type].append(clusterObj.clusterID)
                            RLO.clusterID_count_by_type[cluster_type] += 1
                            # proteins
                            for proteomeID in RLO.proteomeIDs:
                                if proteomeID in clusterObj.proteinIDs_by_proteomeID:
                                    RLO.proteinIDs.extend(clusterObj.proteinIDs_by_proteomeID[proteomeID])
                                    RLO.proteinID_count += clusterObj.proteinID_count_by_proteomeID[proteomeID]
                                    RLO.proteinID_by_type[cluster_type].append(clusterObj.proteinIDs_by_proteomeID[proteomeID])
                                    RLO.proteinID_count_by_type[cluster_type] += clusterObj.proteinID_count_by_proteomeID[proteomeID]
                            # length
                            if self.fasta_parsed == True:
                                for proteinID in RLO.proteinIDs:
                                    RLO.protein_span += sum(clusterObj.protein_length)


    def calculate_rarefraction_data(self, REPETITIONS):
        for rankID in self.proteomeIDs_by_levelID_by_rankID:
            levelIDs = self.levelIDs_by_rankID[rankID]
            for levelID in self.proteomeIDs_by_levelID_by_rankID[rankID]:
                proteomeIDs = self.proteomeIDs_by_levelID_by_rankID[rankID][levelID]
                if len(proteomeIDs) == 1:
                    print "[WARN] - Omitting rarefraction calculation for rank %s at level %s (only contains %s)" % (rankID, levelID, ",".join(proteomeIDs))
                else:
                    print "[STATUS] - Plotting rarefraction curve for rank %s at level %s" % (rankID, levelID)
                    RLO = self.RLO_by_levelID_by_rankID[rankID][levelID]
                    proteomeIDs = self.proteomeIDs_by_levelID_by_rankID[rankID][levelID]
                    for repetition in range(1, REPETITIONS):
                        seen_clusterIDs = set()
                        random_list_of_proteomeIDs = [x for x in proteomeIDs]
                        random.shuffle(random_list_of_proteomeIDs)
                        for idx, proteomeID in enumerate(random_list_of_proteomeIDs):
                            proteome_RLO = self.RLO_by_levelID_by_rankID['proteome'][proteomeID]
                            seen_clusterIDs.update(proteome_RLO.clusterID_by_type['monoton'])
                            seen_clusterIDs.update(proteome_RLO.clusterID_by_type['multiton'])
                            sample_size = idx + 1
                            if not sample_size in RLO.rarefraction_data:
                                RLO.rarefraction_data[sample_size] = []
                            RLO.rarefraction_data[sample_size].append(len(seen_clusterIDs))

    ############################################################################################
    # Writing Output files
    ############################################################################################

    def setup_dirs(self):
        result_path = join(getcwd(), "cb_results")
        self.dirs['main'] = result_path
        print "[STATUS] - Output directories in \n\t%s" % (result_path)
        if not exists(result_path):
            mkdir(result_path)
        for rankID in self.rankIDs:
            rank_path = join(result_path, rankID)
            self.dirs[rankID] = rank_path
            if not exists(rank_path):
                print "\t%s" % (rank_path)
                mkdir(rank_path)

    def output_coverages(self):
        plot_flag = False
        coverages_by_levelID_by_rankID = {}
        for RLO in self.yield_RLOs(ranks=['all'], levels=['all']):
            rankID = RLO.rankID
            levelID = RLO.levelID
            if not rankID in coverages_by_levelID_by_rankID:
                coverages_by_levelID_by_rankID[rankID] = {}
            coverages_by_levelID_by_rankID[rankID][levelID] = [coverage for coverage in sorted(RLO.coverage_in_clusters, reverse=True)]
        for rankID in coverages_by_levelID_by_rankID:
            # PLOT
            coverages_out_png = join(self.dirs[rankID], "%s.coverage.%s" % (rankID, PLOT_FORMAT))
            coverages_by_levelID = coverages_by_levelID_by_rankID[rankID]
            plot_flag = plot_coverage_decay(rankID, coverages_by_levelID, coverages_out_png)
            for levelID in coverages_by_levelID_by_rankID[rankID]:
                # TXT
                if (plot_flag):
                    coverages_out_txt = join(self.dirs[rankID], "%s.%s.coverage.txt" % (rankID, levelID))
                    with open(coverages_out_txt, "w") as fh:
                        fh.write("%s\n" % (coverages_by_levelID_by_rankID[rankID][levelID]))

    def output(self, argument): # for debugging
        if argument == 'categories':
            print "# Ranks : %s" % (self.rankIDs_count)
            print "\t=> %s" % (", ".join(self.rankIDs))
            for rankID in self.rankIDs:
                print "# Levels of rank \"%s\" : %s" % (rankID, self.levelIDs_count_by_rankID[rankID])
                print "\t=> %s" % (", ".join([levelID for levelID in self.levelIDs_by_rankID[rankID]]))
            print "# Proteomes by level by rank"
            for rankID in self.rankIDs:
                print "\t=> %s" % (rankID)
                for level, proteomeIDs in self.proteomeIDs_by_levelID_by_rankID[rankID].items():
                    print "\t %s : %s" % (level, ", ".join(list(proteomeIDs)))
            print "# Levels by rank by proteome"
            for proteomeID in self.proteomeIDs:
                levelIDs = [self.levelIDs_by_rankID_by_proteomeID[proteomeID][rankID] for rankID in self.rankIDs]
                string = []
                for rankID, levelID in zip(self.rankIDs, levelIDs):
                    string.append("\t%s = %s" % (rankID, levelID))
                print "\t - %s \t %s" % (proteomeID, "\t".join(string))
        elif argument == 'ranklevelobjs':
            # all
            #for RLO in self.yield_RLOs(ranks='all', levels='all'):
            #    print RLO
            print "# RLOs"
            for RLO in self.yield_RLOs(ranks=['species'], levels=['all']):
                print RLO
            for RLO in self.yield_RLOs(ranks=['group'], levels=['all']):
                print RLO
            for RLO in self.yield_RLOs(ranks=['genus'], levels=['all']):
                print RLO
        else:
            pass

    def output_counts_by_RLO(self):
        for rankID in self.rankIDs:
            rank_out_f = join(self.dirs[rankID], "%s.counts.txt" % (rankID))
            rank_out_string = ["\t".join(["rankID", "levelID", \
                "cluster_count", "cluster_singleton_count", "cluster_monoton_count", "cluster_multiton_count", \
                "protein_count", "protein_singleton_count", "protein_monoton_count", "protein_multiton_count", \
                "protein_span",  "members_count", "members" \
                ])]
            for levelID, RLO in self.RLO_by_levelID_by_rankID[rankID].items():
                if (self.fasta_parsed):
                    rank_out_string.append("\t".join([str(x) for x in [ \
                        RLO.rankID, RLO.levelID, \
                        RLO.clusterID_count, RLO.clusterID_count_by_type['singleton'], RLO.clusterID_count_by_type['monoton'], RLO.clusterID_count_by_type['multiton'], \
                        RLO.proteinID_count, RLO.proteinID_count_by_type['singleton'], RLO.proteinID_count_by_type['monoton'], RLO.proteinID_count_by_type['multiton'], \
                        RLO.protein_span, RLO.proteomeIDs_count, ",".join(RLO.proteomeIDs) \
                        ]]))
                else:
                    rank_out_string.append("\t".join([str(x) for x in [ \
                        RLO.rankID, RLO.levelID, \
                        RLO.clusterID_count, RLO.clusterID_count_by_type['singleton'], RLO.clusterID_count_by_type['monoton'], RLO.clusterID_count_by_type['multiton'], \
                        RLO.proteinID_count, RLO.proteinID_count_by_type['singleton'], RLO.proteinID_count_by_type['monoton'], RLO.proteinID_count_by_type['multiton'], \
                        "None", RLO.proteomeIDs_count, ",".join(RLO.proteomeIDs) \
                        ]]))
                #level_out_f = join(self.dirs[rankID], "%s.%s.counts.txt" % (rankID, levelID))
                #print levelID
                #print RLO
            with open(rank_out_f, "w") as fh:
                fh.write("\n".join(rank_out_string))

    def output_clusters_by_type(self):
        for rankID in self.rankIDs:
            for levelID, RLO in self.RLO_by_levelID_by_rankID[rankID].items():
                out_cluster_singleton_f = join(self.dirs[rankID], "%s.%s.clusterIDs.singletons.txt" % (rankID, levelID))
                if (RLO.clusterID_by_type['singleton']):
                    with open(out_cluster_singleton_f, "w") as out_cluster_singleton_fh:
                        out_cluster_singleton_fh.write("\n".join(RLO.clusterID_by_type['singleton']))
                out_cluster_monoton_f = join(self.dirs[rankID], "%s.%s.clusterIDs.monotons.txt" % (rankID, levelID))
                if (RLO.clusterID_by_type['monoton']):
                    with open(out_cluster_monoton_f, "w") as out_cluster_monoton_fh:
                        out_cluster_monoton_fh.write("\n".join(RLO.clusterID_by_type['monoton']))
                out_cluster_multiton_f = join(self.dirs[rankID], "%s.%s.clusterIDs.multitons.txt" % (rankID, levelID))
                if (RLO.clusterID_by_type['multiton']):
                    with open(out_cluster_multiton_f, "w") as out_cluster_multiton_fh:
                        out_cluster_multiton_fh.write("\n".join(RLO.clusterID_by_type['multiton']))

    def output_clusters_by_proteome_count(self):
        out_clusterObj_count_f = join(self.dirs['main'], "cluster_count_by_proteome.txt")
        out_clusterObj_count_plot_f = join(self.dirs['main'], "cluster_count_by_proteome.%s" % PLOT_FORMAT)
        out_clusterObj_count_string = ["clusterID\t%s" % "\t".join(self.proteomeIDs)]
        out_clusterObj_count_x = self.proteomeIDs
        out_clusterObj_count_y = []
        out_clusterObj_count_z = []
        for clusterObj in self.yield_clusterObj():
            out_clusterObj_count_y.append(clusterObj.clusterID)
            count_z = [clusterObj.proteinID_count_by_proteomeID.get(proteomeID, 0) for proteomeID in self.proteomeIDs]
            out_clusterObj_count_z.append(count_z)
            out_clusterObj_count_string.append("%s\t%s" % (clusterObj.clusterID, \
                "\t".join([str(clusterObj.proteinID_count_by_proteomeID.get(proteomeID, 0)) for proteomeID in self.proteomeIDs]) \
                ))
        with open(out_clusterObj_count_f, 'w') as out_clusterObj_count_fh:
            out_clusterObj_count_fh.write("\n".join(out_clusterObj_count_string))
        #seaborn_heatmap(out_clusterObj_count_plot_f, out_clusterObj_count_x, out_clusterObj_count_y, out_clusterObj_count_z)
        plot_heatmap(out_clusterObj_count_plot_f, out_clusterObj_count_x, out_clusterObj_count_y, out_clusterObj_count_z)


    def output_rarefraction(self):
        rarefraction_by_levelID_by_rankID = {}
        for RLO in self.yield_RLOs(ranks=['all'], levels=['all']):
            if (RLO.rarefraction_data):
                rankID = RLO.rankID
                levelID = RLO.levelID
                rarefraction_data = RLO.rarefraction_data
                out_f = join(self.dirs[rankID], "%s.%s.rarefraction.txt" % (rankID, levelID))
                if not rankID in rarefraction_by_levelID_by_rankID:
                    rarefraction_by_levelID_by_rankID[rankID] = {}
                rarefraction_by_levelID_by_rankID[rankID][levelID] = rarefraction_data
                print "[STATUS] - Output rarefraction data \t%s" % (out_f)
                with open(out_f, "w") as fh:
                    fh.write("%s\t%s\n" % ("sample_size", "\t".join(["Rep" + str(x) for x in range(1, REPETITIONS)])))
                    for sample_size in rarefraction_data:
                        fh.write("%s\t%s\n" % (sample_size, "\t".join([str(x) for x in rarefraction_data[sample_size]])))
        for rankID in rarefraction_by_levelID_by_rankID:
            rarefraction_plot_f = join(self.dirs[rankID], "%s.rarefraction.%s" % (rankID, PLOT_FORMAT))
            rarefraction_by_levelID = rarefraction_by_levelID_by_rankID[rankID]
            plot_rarefraction_data(rarefraction_by_levelID, rarefraction_plot_f)

############################################################################################
# PLOTTING
############################################################################################

def create_rainbow(ax):
    rainbow = [ax._get_lines.prop_cycler.next()['color']]
    while True:
        nextval = ax._get_lines.prop_cycler.next()['color']
        if nextval not in rainbow:
            rainbow.append(nextval)
        else:
            return rainbow

def next_color(ax):
    rainbow = create_rainbow(ax)
    double_rainbow = collections.deque(rainbow)
    nextval = ax._get_lines.prop_cycler.next()['color']
    double_rainbow.rotate(-1)
    return nextval, itertools.cycle(double_rainbow)

def plot_rarefraction_data(rarefraction_by_levelID, rarefraction_plot_f):
    print "[STATUS] - Plotting rarefraction data \t%s" % (rarefraction_plot_f)
    f, ax = plt.subplots(figsize=(24, 12))
    sns.set_color_codes("pastel")
    order = {}
    #print rarefraction_by_levelID
    #for levelID, rarefraction_data in rarefraction_by_levelID.items():
    #    counter = Counter(rarefraction_data)
    #    order[levelID] = counter[1.0]
    #print order
    max_number_of_samples = 0
    for idx, levelID in enumerate(rarefraction_by_levelID):
        number_of_samples = len(rarefraction_by_levelID[levelID])
        if number_of_samples > max_number_of_samples:
            max_number_of_samples = number_of_samples
        colour = plt.cm.Paired(idx/len(rarefraction_by_levelID))
        #print levelID, colors.rgb2hex(colour)
        x_values = []
        #y_values = []
        y_mins = []
        y_maxs = []
        median_y_values = []
        median_x_values = []

        for x, y_reps in rarefraction_by_levelID[levelID].items():
            #print x, len(y_reps), y_reps
            #x_values.append([x] * len(y_reps))
            #y_values.append(y_reps)
            x_values.append(x)
            y_mins.append(min(y_reps))
            y_maxs.append(max(y_reps))
            median_y_values.append(np.median(y_reps))
            median_x_values.append(x)
        #ax.scatter(x_values, y_values, color=colour, label=levelID, alpha = 0.8)
        x_array = np.array(x_values)
        y_mins_array = np.array(y_mins)
        y_maxs_array = np.array(y_maxs)
        ax.plot(median_x_values, median_y_values, '-', color=colour, label=levelID)
        ax.fill_between(x_array, y_mins_array, y_maxs_array, color=colour, alpha = 0.5)

    ax.set_xlim([0.5, max_number_of_samples + 0.5])
    ax.set_ylabel("Count of clusters")
    ax.set_xlabel("Sampled proteomes")
    ax.legend(ncol=1, numpoints=1, loc="lower right", frameon=True)
    f.savefig(rarefraction_plot_f, format=PLOT_FORMAT)

def cmap_discretize(cmap, N):
    if type(cmap) == str:
        cmap = mat.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in xrange(N+1) ]
    # Return colormap object.
    return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

def plot_heatmap(out_clusterObj_count_plot_f, out_clusterObj_count_x, out_clusterObj_count_y, out_clusterObj_count_z):
    print "[STATUS] - Plotting heatmap \t%s" % (out_clusterObj_count_plot_f)
    # Coordinates
    spacer = 0.02
    dendrogram_height = 0.15
    heatmap_height = heatmap_width = 0.60

    fig_x_min, fig_y_min = 0.1, 0.1
    fig_x_max, fig_y_max = 0.9, 0.9

    rect_row_left = fig_x_min
    rect_row_bottom = fig_y_min
    rect_row_width = dendrogram_height
    rect_row_height = heatmap_height

    rect_col_left = fig_x_min + dendrogram_height + spacer
    rect_col_bottom = fig_y_min + heatmap_height + spacer
    rect_col_width = heatmap_width
    rect_col_height = dendrogram_height + spacer

    rect_heat_left = rect_col_left
    rect_heat_bottom = fig_y_min
    rect_heat_width = heatmap_width
    rect_heat_width = heatmap_height
    rect_row = [rect_row_left, rect_row_bottom, rect_row_width, rect_row_height]
    rect_col = [rect_col_left, rect_col_bottom, rect_col_width, rect_col_height]
    rect_heat = [rect_heat_left,rect_heat_bottom,rect_heat_width,rect_heat_width]

    # setup
    vmin = 0
    vmax = 4
    cmap = cmap_discretize(cm.Greys, N=vmax+1)
    cnorm = mat.colors.Normalize(vmin=vmin, vmax=vmax+1)
    x_labels = out_clusterObj_count_x
    proteome_count = len(x_labels)
    x_values = range(0, proteome_count)
    major_x_ticks = [0.5+x for x in x_values]
    minor_x_ticks = [x for x in x_values]
    FONTSIZE = 32

    fig = pylab.figure(figsize=(48, 48))
    #def median_compare(x, y):
    #    return int(np.mean(y)) - int(np.mean(x))
    #z = np.array([x for x in sorted(out_clusterObj_count_z, cmp=median_compare)])
    z = np.array(out_clusterObj_count_z)

    #print "left started"
    left_Y = sch.linkage(z, method='ward')
    ax_dendrogram_left = fig.add_axes(rect_row, frameon=True)
    #dendrogram_left_Z = sch.dendrogram(left_Y, orientation='left', color_threshold=0.9, show_leaf_counts=True, show_contracted=True, truncate_mode='lastp', distance_sort='ascending')
    dendrogram_left_Z = sch.dendrogram(left_Y, orientation='left', count_sort='descending', color_threshold=0.9, link_color_func=lambda k: '0.5')
    ax_dendrogram_left.set_xticks([])
    ax_dendrogram_left.set_yticks([])
    ax_dendrogram_left.grid(False)
    #print "left done"

    #print "top started"
    top_Y = sch.linkage(z.T, method='ward')
    ax_dendrogram_top = fig.add_axes(rect_col, frameon=True)
    dendrogram_top_Z = sch.dendrogram(top_Y, show_leaf_counts=True, labels=np.array(x_labels), orientation='top', color_threshold=0.9, link_color_func=lambda k: '0.5')
    ax_dendrogram_top.set_yticks([])
    ax_dendrogram_top.set_xticks([])
    ax_dendrogram_top.grid(False)
    #print "top done"

    #reorder leaves
    Z = z[dendrogram_left_Z['leaves'], :]

    axmatrix = fig.add_axes(rect_heat)
    im = axmatrix.matshow(Z, cmap=cmap, aspect='auto', interpolation='nearest', origin='upper', vmin=vmin, vmax=vmax)
    #axmatrix.xaxis.grid(True, which="major")
    axmatrix.xaxis.grid(True)
    axmatrix.set_yticks([])
    #axmatrix.set_xticks(major_x_ticks, minor=False)
    axmatrix.set_xticks([], minor=False)
    axmatrix.set_xticks(minor_x_ticks, minor=True)
    axmatrix.xaxis.grid(False)
    for major_x_tick in major_x_ticks:
        axmatrix.axvline(major_x_tick, linewidth=1, color='1.0')

    axmatrix.xaxis.tick_bottom()
    axmatrix.set_xticklabels(x_labels, minor=True, rotation='vertical', ha='center', fontsize = FONTSIZE)

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.2,0.02,0.5])
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(vmin-0.5, vmax+1.5)
    colorbar = fig.colorbar(mappable, cax=axcolor)
    colorbar.set_ticks(np.linspace(0, vmax+1, vmax+1))
    colorbar.set_ticklabels(["%s" % x if not x == vmax else "%s+" % x for x in range(vmax+1)])
    colorbar.ax.tick_params(labelsize=FONTSIZE)

    fig.savefig(out_clusterObj_count_plot_f)


def plot_coverage_decay(rankID, coverages_by_levelID, coverages_out_png):
    f, ax = plt.subplots(figsize=(24, 12))
    order = {}
    for levelID, coverages in coverages_by_levelID.items():
        counter = Counter(coverages)
        order[levelID] = counter[1.0]

    legend_handles = []
    legend_labels = []
    plot_flag = False
    for levelID_idx, (levelID, count_100) in enumerate(sorted(order.items(), key=operator.itemgetter(1), reverse=True)):
        colour = ''
        label = ''
        number_of_members = dataObj.count_by_levelID_by_rankID[rankID][levelID]
        if number_of_members > 1:
            plot_flag = True

            x_values = []
            y_values = []
            last_x = None
            for coverage_idx, coverage in enumerate(coverages_by_levelID[levelID]):
                if not (y_values):
                    y_values.append(coverage)
                    x_values.append(coverage_idx)
                else:
                    if coverage == y_values[-1]:
                        last_x = coverage_idx
                    else:
                        y_values.append(y_values[-1])
                        x_values.append(last_x)
                        y_values.append(coverage)
                        x_values.append(coverage_idx)
            y_values.append(y_values[-1])
            x_values.append(last_x)
            colour, ax._get_lines.color_cycle = next_color(ax)
            ax.plot(x_values, y_values, '-o', linestyle = '-', linewidth = 4, color = colour, label=levelID, markeredgecolor = 'none')
            label = "%s - %s proteomes" % (levelID, number_of_members)
        if not (colour):
            colour = '0.75'
            label = "%s - %s proteome" % (levelID, number_of_members)
        legend_handles.append(Line2D([0], [0], color=colour, linewidth = 0.5, linestyle="none", marker="o", alpha=1, markersize=10))
        legend_labels.append(label)
    #ax.set_xlim([-0.9, len(x_values)])
    ax.set_ylim([0, 1.1])
    ax.set_ylabel("Fraction of rank members present in clusters")
    ax.set_xlabel("Count of clusters")
    ax.legend(legend_handles, legend_labels, ncol=1, numpoints=1, loc="upper right", frameon=True)
    f.tight_layout()
    if (plot_flag):
        print "[STATUS] - Plotting coverage decay of rankID %s\n\t%s" % (rankID, coverages_out_png)
        f.savefig(coverages_out_png, format=PLOT_FORMAT)
    else:
        print "[WARN] - Omitting coverage decay of rankID %s : too few members in levels " % (rankID)
    return plot_flag

############################################################################################
# CLUSTERS
############################################################################################

class ClusterObj():
    def __init__(self, clusterID, proteinIDs):
        self.clusterID = clusterID
        self.proteinIDs = proteinIDs
        self.proteomeIDs = [x.split(".")[0] for x in proteinIDs]
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
            proteomeID = proteinID.split(".")[0]
            if not proteomeID in proteinIDs_by_proteomeID:
                proteinIDs_by_proteomeID[proteomeID] = set()
            proteinIDs_by_proteomeID[proteomeID].add(proteinID)
        return proteinIDs_by_proteomeID
############################################################################################
# PROTEINS
############################################################################################

class ProteinObj():
    def __init__(self, proteinID, length, proteomeID):
        self.proteinID = proteinID
        self.proteomeID = proteomeID
        self.length = length
        self.clusterID = ''

############################################################################################
# RankLevelObj (PROTEOMES, etc)
############################################################################################

class RankLevelObj():
    def __init__(self, levelID, rank, proteomes):
        self.levelID = levelID
        self.rankID = rank
        self.proteomeIDs = proteomes
        self.proteomeIDs_count = len(proteomes)
        self.file = None

        self.proteinIDs = []
        self.proteinID_count = 0
        self.proteinID_by_type = {'singleton' : [], 'monoton' : [], 'multiton' : []}
        self.proteinID_count_by_type = {'singleton' : 0, 'monoton' : 0, 'multiton' : 0}
        self.protein_span = 0
        self.clusterIDs = []
        self.clusterID_count = 0
        self.clusterID_by_type = {'singleton' : [], 'monoton' : [], 'multiton' : []}
        self.clusterID_count_by_type = {'singleton' : 0, 'monoton' : 0, 'multiton' : 0}

        self.coverage_in_clusters = []
        self.rarefraction_data = {} # repetition : number of clusters

    def __str__(self):
        string = ''
        string += "levelID : %s\n" % self.levelID
        string += "rankID : %s\n" % self.rankID
        string += "proteinID_count : %s\n" % self.proteinID_count
        string += "proteinID_count_by_type : %s\n" % self.proteinID_count_by_type
        string += "clusterID_count : %s\n" % self.clusterID_count
        string += "clusterID_count_by_type : %s\n" % (str(self.clusterID_count_by_type))
        string += "coverage : %s\n" % (str(Counter(self.coverage_in_clusters)))
        return string

    def check(self):
        count_monoton_and_multiton_clusters = self.clusterID_count_by_type.get('monoton', 0) + self.clusterID_count_by_type.get('multiton', 0)
        count_coverages = len(self.coverage_in_clusters)
        if count_monoton_and_multiton_clusters != count_coverages:
            return False
        elif self.proteinID_count != len(set(self.proteinIDs)):
            return False
        elif self.clusterID_count != len(set(self.clusterIDs)):
            return False
        else:
            return True

############################################################################################

if __name__ == "__main__":
    # david aanensen, github?
    REPETITIONS = 30 + 1
    PLOT_FORMAT = 'png'
    species_ids_f, species_classification_f, groups_f = '','',''
    try:
        species_ids_f = sys.argv[1]
        species_classification_f = sys.argv[2]
        groups_f = sys.argv[3]
    except:
        sys.exit("./clusterbuster.py SPECIESID_FILE SPECIESCLASSIFICATION_FILE GROUPS_FILE")

    dataObj = DataObj()
    # Get all info from species_classification_f
    dataObj.parse_species_classification(species_classification_f)
    dataObj.create_RLOs()
    dataObj.parse_species_ids(species_ids_f)
    dataObj.setup_dirs()
    dataObj.parse_fasta()
    #dataObj.output("categories") # debug
    dataObj.parse_clusters(groups_f)
    #dataObj.output("ranklevelobjs")
    dataObj.output_coverages()
    dataObj.calculate_rarefraction_data(REPETITIONS)
    dataObj.output_rarefraction()
    dataObj.output_counts_by_RLO()
    dataObj.output_clusters_by_type()
    dataObj.output_clusters_by_proteome_count()
    #for RLO in dataObj.yield_RLOs(ranks=['all'], levels=['all']):
    #    print RLO.rankID, RLO.levelID, RLO.check()
