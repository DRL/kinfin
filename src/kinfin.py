#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
import os
import numpy as np
import itertools
import random

import matplotlib as mat
import matplotlib.cm as cm
mat.use('agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

'''
RULES:
- All internal keys are prefix_id and cluster_id

To Do:
- remove formating for fractions for clusters by species
- add counts for clusters by species


- classification_f:
categories : columns
elements : elements of categories
ulements : unique elements

Data:
A) plot 1:1 shared by X taxa:
    - for each inflation value
        - for each category (proteome, species, genus, ...)
            - for each combination of X taxa within category (ulements)
                    - get number of 1:1's limited to those taxa

B) plot Rarefraction curve:
    - for each inflation value
        - for each category (proteome, species, genus, ...)
            - get rarefraction data for N repetitions

C) plot Cluster/Protein Count by type (VENN)
    1. cluster/protein count by type for each category
        for each inflation value
            - for each PROTEOME
                -  for each category
                    - get count of clusters/proteins that are:
                        - singletons
                        - monotons
                        - multitons with any other ulement in category
                        => CP1207 :
                            -singleton cluster/protein count
                            -monoton cluster/protein count
                            -multiton cluster/protein count w/ CMPRLs
                            -multiton cluster/protein count w/ notCMPRLs
                            -multiton cluster/protein count w/ CBOT
                            -multiton cluster/protein count w/ CBOT + CMPRLs
                            -multiton cluster/protein count w/ CBOT + CMPRLs

D) For each inflation value:
    - for each PROTEOME
        -  for each category
            get count proteins/clusters
                - singletons
                - monotons (same species)

            pan = core (in all strains) + accessory (not in all strains)
            accessory = singletons + monotons + multi-within + multi-with-others = not core


'''


#############################################
################## CLASSES ##################
#############################################
class OldCladeObj():
    def __init__(self, name):
        self.name = name
        self.cluster_multispecies_count = 0
        self.protein_multispecies_count = 0
        self.cluster_monospecies_count = 0
        self.protein_monospecies_count = 0
        self.cluster_singleton_count = 0
        self.protein_singleton_count = 0

# DATA
class OldDataObj():
    def __init__(self, inflation_value):
        self.order_of_clusters = [] # IDs
        self.order_of_species_ids = [] # IDs
        self.set_of_species_names = set() # names
        self.set_of_nematode_names = set()
        self.set_of_nematode_ids = set()
        self.inflation_value = inflation_value
        self.species_by_id = {} # by ID
        self.clusters_by_id = {} # by ID

        self.one2one_clusters = {} # by percentage

        self.protein_count = 0
        self.cluster_count = 0

        self.cluster_multispecies_count = 0
        self.protein_multispecies_count = 0

        self.cluster_monospecies_count = 0
        self.protein_monospecies_count = 0

        self.protein_singleton_count = 0

    # CALCULATIONS
    def calc_median_cluster_size(self, cluster_type, by):
        medians = []
        for clusterObj in self.yield_clusterObjs():
            if cluster_type == clusterObj.type:
                if by == "id":
                    medians.append(clusterObj.proteins_by_species_id_median)
                elif by == "name":
                    medians.append(clusterObj.proteins_by_species_name_median)
        if medians:
            return np.median(medians)
        else:
            return 0

    def calc_mean_cluster_size(self, cluster_type, by):
        means = []
        for clusterObj in self.yield_clusterObjs():
            if cluster_type == clusterObj.type:
                if by == "id":
                    means.append(clusterObj.proteins_by_species_id_mean)
                elif by == "name":
                    means.append(clusterObj.proteins_by_species_name_mean)

        if (means):
            return np.mean(means)
        else:
            return 0


    # YIELD
    def yield_speciesObjs(self):
        for order in self.order_of_species_ids:
            yield self.species_by_id[order]

    def yield_clusterObjs(self):
        for order in self.order_of_clusters:
            yield self.clusters_by_id[order]

    # output
    def output_rarefraction_data_by_clade(self, population, repetitions):
        if population == "nematode":
            nematode_clades = set([x for x in CLASS.values() if not x == "Outgroup"])
            for clade in nematode_clades:
                data_all = {}
                #print "clade: %s" % (clade)
                for repetition in range(1, repetitions):
                    #print "repetition: %s" % (repetition)
                    seen_cluster_ids = set()
                    sample_size = 0

                    random_list_of_nematode_ids_in_clade = list(x for x in self.set_of_nematode_ids if self.species_by_id[x].clade == clade)
                    random.shuffle(random_list_of_nematode_ids_in_clade)
                    for nematode_id in random_list_of_nematode_ids_in_clade:
                        sample_size += 1
                        speciesObj = self.species_by_id[nematode_id]
                        seen_cluster_ids.update(speciesObj.cluster_ids)
                        if not sample_size in data_all:
                            data_all[sample_size] = []
                        #print "sample_size: %s" % (sample_size), "cluster_count: ", len(seen_cluster_ids)
                        data_all[sample_size].append(len(seen_cluster_ids))

                rarefraction_stats_f = self.inflation_value + "." + clade + ".rarefraction_stats.txt"
                with open(rarefraction_stats_f, "w") as fh:
                    fh.write("%s\t%s\t%s\t%s\t%s\n" % (\
                        "sample_size", \
                        "clusters_data_all_mean", \
                        "clusters_data_all_max", \
                        "clusters_data_all_min", \
                        "\t".join(["Rep" + str(x) for x in range(1, repetitions)]) \
                        ))

                    for sample_size in sorted(data_all):
                        clusters_data_all = data_all[sample_size]
                        clusters_data_all_mean = np.mean(clusters_data_all)
                        clusters_data_all_max = np.max(clusters_data_all)
                        clusters_data_all_min = np.min(clusters_data_all)

                        fh.write("%s\t%s\t%s\t%s\t%s\n" % (\
                            sample_size, \
                            ("%.2f" % (clusters_data_all_mean)), \
                            clusters_data_all_max, \
                            clusters_data_all_min, \
                            "\t".join([str(x) for x in clusters_data_all]) \
                            ))

    def output_rarefraction_data(self, population, repetitions):
        if population == "nematode":
            data_all = {}

            count_of_nematode_ids = len(self.set_of_nematode_ids)
            #print "number of nematode ids: %s" % (count_of_nematode_ids)
            for repetition in range(1, repetitions):
                #print "repetition: %s" % (repetition)
                # bootstrapping

                seen_cluster_ids = set()
                sample_size = 0

                random_list_of_nematode_ids = list(self.set_of_nematode_ids)
                random.shuffle(random_list_of_nematode_ids)

                for nematode_id in random_list_of_nematode_ids:
                    sample_size += 1
                    speciesObj = self.species_by_id[nematode_id]
                    seen_cluster_ids.update(speciesObj.cluster_ids)
                    if not sample_size in data_all:
                        data_all[sample_size] = []
                    #print "sample_size: %s" % (sample_size), "cluster_count: ", len(seen_cluster_ids)
                    data_all[sample_size].append(len(seen_cluster_ids))

            rarefraction_stats_f = self.inflation_value + ".rarefraction_stats.txt"
            with open(rarefraction_stats_f, "w") as fh:
                fh.write("%s\t%s\t%s\t%s\t%s\n" % (\
                    "sample_size", \
                    "clusters_data_all_mean", \
                    "clusters_data_all_max", \
                    "clusters_data_all_min", \
                    "\t".join(["Rep" + str(x) for x in range(1, repetitions)]) \
                    ))

                for sample_size in sorted(data_all):
                    clusters_data_all = data_all[sample_size]
                    clusters_data_all_mean = np.mean(clusters_data_all)
                    clusters_data_all_max = np.max(clusters_data_all)
                    clusters_data_all_min = np.min(clusters_data_all)

                    fh.write("%s\t%s\t%s\t%s\t%s\n" % (\
                        sample_size, \
                        ("%.2f" % (clusters_data_all_mean)), \
                        clusters_data_all_max, \
                        clusters_data_all_min, \
                        "\t".join([str(x) for x in clusters_data_all]) \
                        ))
        else:
            sys.exit("[ERROR7]")


    def output_one2one_clusters_by_id(self, percentages):
        for percentage in percentages:
            output = ""
            for clusterObj in self.yield_clusterObjs():
                if clusterObj.nematode_id_fraction >= percentage:
                    good = 1
                    for species_id, protein_count in clusterObj.proteins_by_species_id.items():
                        if protein_count > 1:
                            good = 0
                    if (good):
                        output += "%s\t%s\t%s\n" % ( \
                            clusterObj.id, \
                            clusterObj.nematode_id_fraction, \
                            ",".join(clusterObj.proteins))
            if (output):
                with open(str(self.inflation_value) + "_one2one_" + str(percentage) + ".species_id.txt", "w") as fh:
                    fh.write("%s\t%s\t%s\n" % ("cluster_id", "nematode_id_fraction", "proteins"))
                    fh.write(output)



    def output_one2one_clusters_by_name(self, percentages):
        for percentage in percentages:
            output = ""
            for clusterObj in self.yield_clusterObjs():
                if clusterObj.nematode_name_fraction >= percentage:
                    print clusterObj.__dict__
                    good = 1
                    for species_name, protein_count in clusterObj.proteins_by_species_name.items():
                        print species_name, protein_count
                        print assemblies
                        if protein_count > assemblies[species_name]: # based on how many assemblies are present for a species
                            good = 0
                    if (good):
                        output += "%s\t%s\t%s\n" % ( \
                            clusterObj.id, \
                            clusterObj.nematode_name_fraction, \
                            ",".join(clusterObj.proteins))
            if (output):
                with open(str(self.inflation_value) + "_one2one_" + str(percentage) + ".species_name.txt", "w") as fh:
                    fh.write("%s\t%s\t%s\n" % ("cluster_id", "nematode_id_fraction", "proteins"))
                    fh.write(output)

    # ADDING SPECIES
    def add_speciesObj(self, speciesObj):
        self.order_of_species_ids.append(speciesObj.id)
        self.set_of_species_names.add(speciesObj.name)
        self.species_by_id[speciesObj.id] = speciesObj
        if not speciesObj.clade == 'Outgroup':
            self.set_of_nematode_ids.add(speciesObj.id)
            self.set_of_nematode_names.add(speciesObj.name)

    # ADDING CLUSTERS
    def add_clusterObj(self, clusterObj):

        #   Update ClusterObj
        ##  species fraction of cluster
        clusterObj.species_id_fraction = clusterObj.species_id_count / len(self.order_of_species_ids)
        clusterObj.species_name_fraction = clusterObj.species_name_count / len(self.set_of_species_names)
        # add fraction
        nematode_ids_in_cluster = set()
        nematode_names_in_cluster = set()
        for sp_id in clusterObj.species_ids:
            if not CLASS[sp_id] == 'Outgroup':
                nematode_ids_in_cluster.add(sp_id)
                nematode_names_in_cluster.add(self.species_by_id[sp_id].name)

        clusterObj.nematode_id_fraction = len(nematode_ids_in_cluster) / len(self.set_of_nematode_ids)
        clusterObj.nematode_name_fraction = len(nematode_names_in_cluster) / len(self.set_of_nematode_names)


        # Add clusterObj
        self.order_of_clusters.append(clusterObj.id)
        self.clusters_by_id[clusterObj.id] = clusterObj

        self.cluster_count += 1
        self.protein_count += clusterObj.protein_count

        # set cluster types and sum for DataObj
        if clusterObj.species_name_count > 1:
            self.cluster_multispecies_count += 1
            self.protein_multispecies_count += clusterObj.protein_count
            clusterObj.type = 'multispecies'
        elif clusterObj.species_name_count == 1:
            if clusterObj.protein_count > 1:
                self.cluster_monospecies_count += 1
                self.protein_monospecies_count += clusterObj.protein_count
                clusterObj.type = 'monospecies'
            elif clusterObj.protein_count == 1:
                self.protein_singleton_count += 1
                clusterObj.type = 'singleton'
            else:
                sys.exit("[ERROR 6]")
        else:
            sys.exit("[ERROR 7]")

        #   Update SpeciesObj
        for species_id, protein_count in clusterObj.proteins_by_species_id.items():

            ##  Cluster-membership for each species
            self.species_by_id[species_id].cluster_ids.add(clusterObj.id)
            self.species_by_id[species_id].cluster_count += 1

            ##  Protein counts for each species
            self.species_by_id[species_id].protein_count += protein_count

            ## Protein counts for multi-species-clusters (based on species_name)
            if clusterObj.type == 'multispecies':
                self.species_by_id[species_id].cluster_multispecies_count += 1
                self.species_by_id[species_id].cluster_multispecies_ids.add(clusterObj.id)
                self.species_by_id[species_id].protein_multispecies_count += protein_count
            elif clusterObj.type == 'monospecies':
                self.species_by_id[species_id].cluster_monospecies_count += 1
                self.species_by_id[species_id].cluster_monospecies_ids.add(clusterObj.id)
                self.species_by_id[species_id].protein_monospecies_count += protein_count
            elif clusterObj.type == 'singleton':
                self.species_by_id[species_id].protein_singleton_count += 1
            else:
                sys.exit("[ERROR] - 5")

    # PARSING SPECIES IDS
    def parse_species_ids(self, species_ids_f):
        with open(species_ids_f) as fh:
            for line in fh:
                number, species_f = line.rstrip("\n").split(": ")
                speciesObj = SpeciesObj(species_f)
                self.add_speciesObj(speciesObj)

    # PARSING GROUPS
    def parse_groups(self, groups_f):
        i = 0
        with open(groups_f) as fh:
            for line in fh:
                i += 1
                if (i % 1000) == 0:
                    sys.stdout.write('\r')
                    print i,
                    sys.stdout.flush()
                cluster_id, protein_string = line.rstrip("\n").split(": ")
                clusterObj = ClusterObj(cluster_id, protein_string.split(), self.species_by_id)
                self.add_clusterObj(clusterObj)
        print


# CLUSTERS
class OldClusterObj():
    def __init__(self, cluster_id, proteins, species_by_id):
        self.id = cluster_id
        self.proteins = set(proteins)
        self.protein_count = len(proteins)

        self.proteins_by_species_id = self.calc_proteins_by_species_id(proteins)
        self.species_ids = set(self.proteins_by_species_id.keys())
        self.species_id_count = len(self.species_ids)
        self.proteins_by_species_id_median = np.median([v for v in self.proteins_by_species_id.values()])
        self.proteins_by_species_id_mean = np.mean([v for v in self.proteins_by_species_id.values()])

        self.proteins_by_species_name = self.calc_proteins_by_species_name(proteins, species_by_id)
        self.species_names = set(self.proteins_by_species_name.keys())
        self.species_name_count = len(self.species_names)
        self.proteins_by_species_name_median = np.median([v for v in self.proteins_by_species_name.values()])
        self.proteins_by_species_name_mean = np.mean([v for v in self.proteins_by_species_name.values()])

        # changed later
        self.species_id_fraction = 0.0
        self.species_name_fraction = 0.0
        self.nematode_id_fraction = 0.0
        self.nematode_name_fraction = 0.0
        self.type = ''

    def get_proteins_of_species(self, species_id):
        proteins_of_species_id = []
        for protein in self.proteins:
            if species_id == protein.split(".")[0]:
                proteins_of_species_id.append(protein)
        return proteins_of_species_id

    def calc_proteins_by_species_id(self, proteins):
        '''
        - only stores dict for species ids in that cluster
        '''
        temp = {}
        for protein in proteins:
            species_id = protein.split(".")[0]
            temp[species_id] = temp.get(species_id, 0) + 1
        return temp

    def calc_proteins_by_species_name(self, proteins, species_by_id):
        '''
        - only stores dict for species names in that cluster
        '''
        temp = {}
        for protein in proteins:
            species_id = protein.split(".")[0]
            species_name = species_by_id[species_id].name
            temp[species_name] = temp.get(species_name, 0) + 1
        return temp

#############################################
################## FUNCTIONS ################
#############################################

def output(results):

    print "Writing results..."

    group_stats_fh = open("group_stats.txt", "w")
    #group_stats_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( \
    group_stats_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( \
            "inflation_value", \
            "total_species_id_count", \
            "total_species_name_count", \
            "cluster_count", \
            "protein_count", \
            "cluster_singleton_count", \
            "protein_singleton_count", \
            "cluster_singleton_fraction", \
            "protein_singleton_fraction", \
            "cluster_monospecies_count", \
            "protein_monospecies_count", \
            "cluster_monospecies_fraction", \
            "protein_monospecies_fraction", \
            "cluster_multispecies_count", \
            "protein_multispecies_count", \
            "cluster_multispecies_fraction", \
            "protein_multispecies_fraction", \
            "cluster_multispecies_by_name_mean", \
            "cluster_multispecies_by_name_median", \
            "cluster_multispecies_by_id_mean", \
            "cluster_multispecies_by_id_median"))

    species_stats_fh = open("species_stats.txt", "w")
    species_stats_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( \
        "inflation_value", \
        "species_id", \
        "species_name", \
        "species_source", \
        "species_type", \
        "species_class", \
        "species_cluster_count", \
        "species_protein_count", \
        "species_cluster_singleton_fraction", \
        "species_protein_singleton_fraction", \
        "species_cluster_monospecies_fraction", \
        "species_protein_monospecies_fraction", \
        "species_cluster_multispecies_fraction", \
        "species_protein_multispecies_fraction"))

    cluster_stats_fh = open("cluster_stats.txt", "w")
    cluster_stats_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        "inflation_value", \
        "cluster_id", \
        "cluster_protein_count", \
        "cluster_species_id_count", \
        "cluster_species_name_count", \
        "cluster_median_cluster_size", \
        "cluster_mean_cluster_size", \
        "cluster_species_id_fraction", \
        "cluster_species_name_fraction", \
        "cluster_proteins"))

    inflation_values = sorted(results.keys())
    species_list = results[inflation_values[0]].order_of_species_ids
    cluster_binary_fh = open("cluster_binary_stats.txt", "w")
    cluster_binary_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" % ( \
        "inflation_value", \
        "cluster_id", \
        "\t".join(species_list),\
        "protein_count",\
        "cluster_median_cluster_size", \
        "cluster_mean_cluster_size", \
        "cluster_species_id_fraction", \
        "cluster_species_name_fraction"))


    print inflation_values

    for inflation_value in inflation_values:
        print inflation_value
        dataObj = results[inflation_value]

        # sampling curve
        print "Calculation of collector's curve data"
        dataObj.output_rarefraction_data("nematode", 25)
        dataObj.output_rarefraction_data_by_clade("nematode", 25)
        total_species_id_count = len(dataObj.order_of_species_ids)
        total_species_name_count = len(dataObj.set_of_species_names)

        protein_count = dataObj.protein_count

        cluster_count = dataObj.cluster_count

        cluster_multispecies_count = dataObj.cluster_multispecies_count
        protein_multispecies_count = dataObj.protein_multispecies_count
        cluster_multispecies_fraction = '{:.2}'.format(dataObj.cluster_multispecies_count/cluster_count)
        protein_multispecies_fraction = '{:.2}'.format(dataObj.protein_multispecies_count/protein_count)

        cluster_monospecies_count = dataObj.cluster_monospecies_count
        protein_monospecies_count = dataObj.protein_monospecies_count
        cluster_monospecies_fraction = '{:.2}'.format(dataObj.cluster_monospecies_count/cluster_count)
        protein_monospecies_fraction = '{:.2}'.format(dataObj.protein_monospecies_count/protein_count)

        cluster_singleton_count = dataObj.protein_singleton_count
        protein_singleton_count = dataObj.protein_singleton_count
        cluster_singleton_fraction = '{:.2}'.format(dataObj.protein_singleton_count/cluster_count)
        protein_singleton_fraction = '{:.2}'.format(dataObj.protein_singleton_count/protein_count)

        cluster_multispecies_by_name_mean = dataObj.calc_mean_cluster_size("multispecies", "name")
        #cluster_monospecies_by_name_mean = dataObj.calc_mean_cluster_size("monospecies", "name")
        #cluster_singleton_by_name_mean = dataObj.calc_mean_cluster_size("singletons", "name")

        cluster_multispecies_by_id_mean = dataObj.calc_mean_cluster_size("multispecies", "id")
        #cluster_monospecies_by_id_mean = dataObj.calc_mean_cluster_size("monospecies", "id")
        #cluster_singleton_by_id_mean = dataObj.calc_mean_cluster_size("singletons", "id")

        cluster_multispecies_by_name_median = dataObj.calc_median_cluster_size("multispecies", "name")
        #cluster_monospecies_by_name_median = dataObj.calc_median_cluster_size("monospecies", "name")
        #cluster_singleton_by_name_median = dataObj.calc_median_cluster_size("singletons", "name")

        cluster_multispecies_by_id_median = dataObj.calc_median_cluster_size("multispecies", "id")
        #cluster_monospecies_by_id_median = dataObj.calc_median_cluster_size("monospecies", "id")
        #cluster_singleton_by_id_median = dataObj.calc_median_cluster_size("singletons", "id")
        print "Writing group stats"
        group_stats_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( \
            inflation_value, \
            total_species_id_count, \
            total_species_name_count, \
            cluster_count, \
            protein_count, \
            cluster_singleton_count, \
            protein_singleton_count, \
            cluster_singleton_fraction, \
            protein_singleton_fraction, \
            cluster_monospecies_count, \
            protein_monospecies_count, \
            cluster_monospecies_fraction, \
            protein_monospecies_fraction, \
            cluster_multispecies_count, \
            protein_multispecies_count, \
            cluster_multispecies_fraction, \
            protein_multispecies_fraction, \
            cluster_multispecies_by_name_mean, \
            cluster_multispecies_by_name_median, \
            cluster_multispecies_by_id_mean, \
            cluster_multispecies_by_id_median))

        # Species by inflation_value
        print "Writing species stats"
        for speciesObj in dataObj.yield_speciesObjs():

            species_id = speciesObj.id
            species_name = speciesObj.name
            species_source = speciesObj.source
            species_type = speciesObj.type
            species_class = speciesObj.clade
            species_cluster_count = speciesObj.cluster_count
            species_protein_count = speciesObj.protein_count

            if species_cluster_count:
                species_cluster_singleton_fraction = '{:.2}'.format(speciesObj.protein_singleton_count/species_cluster_count)
                species_protein_singleton_fraction = '{:.2}'.format(speciesObj.protein_singleton_count/species_protein_count)

                species_cluster_monospecies_fraction = '{:.2}'.format(speciesObj.cluster_monospecies_count/species_cluster_count)
                species_protein_monospecies_fraction = '{:.2}'.format(speciesObj.protein_monospecies_count/species_protein_count)

                species_cluster_multispecies_fraction = '{:.2}'.format(speciesObj.cluster_multispecies_count/species_cluster_count)
                species_protein_multispecies_fraction = '{:.2}'.format(speciesObj.protein_multispecies_count/species_protein_count)
            else:
                species_cluster_singleton_fraction = '{:.2}'.format(0.0)
                species_protein_singleton_fraction = '{:.2}'.format(0.0)
                species_cluster_monospecies_fraction = '{:.2}'.format(0.0)
                species_protein_monospecies_fraction = '{:.2}'.format(0.0)
                species_cluster_multispecies_fraction = '{:.2}'.format(0.0)
                species_protein_multispecies_fraction = '{:.2}'.format(0.0)

            species_stats_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                inflation_value, \
                species_id, \
                species_name, \
                species_source, \
                species_type, \
                species_class, \
                species_cluster_count, \
                species_protein_count, \
                species_cluster_singleton_fraction, \
                species_protein_singleton_fraction, \
                species_cluster_monospecies_fraction, \
                species_protein_monospecies_fraction, \
                species_cluster_multispecies_fraction, \
                species_protein_multispecies_fraction))

        # Clusters by inflation_value
        print "Writing cluster/binary stats"
        for clusterObj in dataObj.yield_clusterObjs():

            cluster_id = clusterObj.id
            cluster_protein_count = clusterObj.protein_count
            cluster_species_id_count = clusterObj.species_id_count
            cluster_species_name_count = clusterObj.species_name_count
            cluster_median_cluster_size = clusterObj.proteins_by_species_name_median
            cluster_mean_cluster_size = clusterObj.proteins_by_species_name_mean
            cluster_species_id_fraction = clusterObj.species_id_fraction
            cluster_species_name_fraction = clusterObj.species_name_fraction

            cluster_proteins = ",".join(clusterObj.proteins)
            cluster_stats_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                inflation_value, \
                cluster_id, \
                cluster_protein_count, \
                cluster_species_id_count, \
                cluster_species_name_count, \
                cluster_median_cluster_size, \
                cluster_mean_cluster_size, \
                cluster_species_id_fraction, \
                cluster_species_name_fraction, \
                cluster_proteins))

            binary_string = ["1" if (sp in clusterObj.species_ids) else "0" for sp in species_list]
            cluster_binary_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" % ( \
                inflation_value, \
                cluster_id, \
                "\t".join(binary_string),\
                cluster_protein_count,\
                cluster_median_cluster_size, \
                cluster_mean_cluster_size, \
                cluster_species_id_fraction, \
                cluster_species_name_fraction))

        print "Writing one2one stats"
        dataObj.output_one2one_clusters_by_id(PERCENTAGES)
        dataObj.output_one2one_clusters_by_name(PERCENTAGES)


    species_stats_fh.close()
    group_stats_fh.close()
    cluster_stats_fh.close()
    cluster_binary_fh.close()

def parse_classification(species_classification_f):
    clade = {}
    assemblies = {}
    with open(species_classification_f) as fh:
        for l in fh:
            line = l.rstrip("\n").split()
            clade[line[0]] = line[1]
            assemblies[line[2]] = assemblies.get(line[2], 0) + 1
    return clade, assemblies

def generate_dataObjs(species_ids_f, groups_fs):
    results = {}
    for groups_f in groups_fs:
        if os.path.isfile(groups_f) and groups_f.endswith(".txt"):
            print "Parsing %s" % (groups_f)
            inflation_value = os.path.basename(groups_f).lstrip("OrthologousGroups_").rstrip(".txt")
            dataObj = DataObj(inflation_value)
            dataObj.parse_species_ids(species_ids_f)
            dataObj.parse_groups(groups_f)
            results[inflation_value] = dataObj
    return results

################ NEW

# Data
class DataObj():
    def __init__(self):
        self.proteomeObjs = {} # by ID
        self.proteome_count = 0
        self.proteome_order = []

        self.proteins = {} # proteinID : clusterID
        self.clusterObjs = {} # by IV
        self.cluster_count = {} # by IV
        self.cluster_order = [] # by IV

        '''
        self.category_keys: is read from header of species_classification_f
        self.category_labels: is read from rows of species_classification_f and unique'd
        self.category_labels_combinations: is build from taking all combinations of self.category_labels

        => These are used to get counts (cluster, proteins) and cluster names of all subsets for each category (venn diagram)
        '''

        self.category_keys = []
        self.category_labels = {}
        self.category_labels_combinations = set()

        self.inflation_values = []
        self.dirs = {}

    def add_proteomeObjs(self, species_ids_f):
        with open(species_ids_f) as fh:
            for l in fh:
                number, species_f = l.rstrip("\n").split(": ")
                proteomeObj = ProteomeObj(number, species_f)
                self.proteomeObjs[proteomeObj.id] = proteomeObj
                self.proteome_count += 1
                self.proteome_order.append(proteomeObj.id)

    def add_ranks_to_proteomeObjs(self, categoryObj):
        for proteomeObj_id, ranks in categoryObj.ranks_by_proteome.items():
            self.proteomeObjs[proteomeObj_id].ranks = ranks


    def setup_dirs(self):
        result_path = os.path.join(os.getcwd(), "cb_results")
        self.dirs['results'] = result_path
        print "[STATUS] - Creating directories \n\t%s" % (result_path)
        if not os.path.exists(result_path):
            os.mkdir(result_path)
        for label in self.label_keys:
            label_path = os.path.join(result_path, label)
            self.dirs[label] = label_path
            if not os.path.exists(label_path):
                print "\t%s" % (label_path)
                os.mkdir(label_path)

    def add_orthologous_groups(self, groups_fs):
        for groups_f in groups_fs:
            if os.path.isfile(groups_f) and groups_f.endswith(".txt"):
                inflation_value = os.path.basename(groups_f).lstrip("OrthologousGroups_").rstrip(".txt")
                print "[STATUS] - Parsing %s : inflation value = %s" % (groups_f, inflation_value)
                self.parse_group_f(inflation_value, groups_f)
            else:
                sys.exit("[ERROR] - %s is not a file" % (groups_f))

    def parse_group_f(self, inflation_value, groups_f):
        with open(groups_f) as fh:
            for line in fh:
                cluster_id, protein_string = line.rstrip("\n").split(": ")
                clusterObj = ClusterObj(cluster_id, protein_string.split(), inflation_value)
                self.add_clusterObj(clusterObj)

    def add_clusterObj(self, clusterObj):
        '''
         ClusterObjs are created
          - protein count
          - for each category
                - set of members
                - count of members
                - type of cluster : singleton/monoton/multiton
                - if not singleton (singleton can't be core)
                     - get percentage of total of group (store where?)


        '''


        for idx, classification in enumerate(self.categories):
            flavours_in_cluster = [self.proteomeObjs[proteome_id].categories[idx] for proteome_id in clusterObj.proteomes]
            clusterObj.protein_count_by_class[classification] = {}
            for flavour in flavours_in_cluster:
                clusterObj.protein_count_by_class[classification][flavour] = clusterObj.protein_count_by_class[classification].get(flavour, 0) + 1
            unique_flavours_in_cluster = set(flavours_in_cluster)
            # cluster type by class
            if clusterObj.protein_count == 1:
                clusterObj.type_by_class[classification] = "singleton"
            else:
                if len(unique_flavours_in_cluster) == 1:
                    clusterObj.type_by_class[classification] = "mono"
                else:
                    clusterObj.type_by_class[classification] = "multi"
        print clusterObj.__dict__

        #clusterObj.proteomes
        if not clusterObj.inflation_value in self.inflation_values:
            self.inflation_values.append(clusterObj.inflation_value)
            self.clusterObjs[clusterObj.inflation_value] = {}
            self.cluster_count[clusterObj.inflation_value] = 0
            self.cluster_order[clusterObj.inflation_value] = []
        self.clusterObjs[clusterObj.inflation_value][clusterObj.id] = clusterObj
        self.cluster_count[clusterObj.inflation_value] += 1
        self.cluster_order[clusterObj.inflation_value].append(clusterObj.id)

    def yield_clusterObjs(self):
        '''
        yields all clusterObjs across all inflation values
        '''
        for inflation_value in self.inflation_values:
            for cluster_id in self.cluster_order[inflation_value]:
                yield self.clusterObjs[inflation_value][cluster_id]

    def yield_proteomeObjs(self):
        for proteome_id in self.proteome_order:
            yield self.proteomeObjs[proteome_id]

    def update_clusterObjs(self):
        for classification in self.categories:
            for clusterObj in self.yield_clusterObjs():
                pass



# CLUSTERS
class ClusterObj():
    def __init__(self, cluster_id, proteins, inflation_value):
        self.id = cluster_id
        self.inflation_value = inflation_value
        self.proteins = set(proteins)
        self.protein_count = len(proteins)
        self.proteomes = [x.split(".")[0] for x in proteins]

        '''
        member_count_at_rank['species'] = ('X' : 2, 'O' : 3)
        member_count_at_rank['proteomes'] = ('A' : 2, 'B' : 1, 'C' : 1,'D' : 1, 'E' :1)
        members_at_rank['species'] = set(member_count_at_rank['species'].keys())        # ('X', 'O')
        members_at_rank['proteome'] = set(member_count_at_rank['proteomes']).keys())    # ('A', 'B', 'C','D', 'E')
        member_perc_at_rank['species'] = ('X' : 1.0, 'O' : 1.0)
        member_perc_at_rank['proteomes'] = ('A' : 1.0, 'B' : 1.0, 'C' : 1.0,'D' : 1.0, 'E' :1.0)
        cluster_type_at_rank['species'] = (multiton)
        cluster_type_at_rank['proteome'] = (multiton)
        '''

        { : for proteome in self.proteomes}
        self.members_at_rank = {}
        self.member_count_at_rank = {}
        self.member_perc_at_rank = {}
        self.cluster_type_at_rank = {}

        self.category_of_members = {}

        self.type_by_class = {}
        self.protein_count_by_class = {}
        self.fraction_by_class = {}
        self.mean_by_class = {}
        self.median_by_class = {}

# Proteomes
class ProteomeObj():
    def __init__(self, number, species_f):
        self.file = species_f
        self.fields = species_f.split(".")
        self.id = self.fields[0]
        self.idx = int(number)
        self.ranks = set()

        ## changed later
        #self.cluster_ids = set() # cluster_ids of clusters the species is a member of
        #self.cluster_count = 0
        #self.protein_count = 0
        #self.cluster_monospecies_ids = set()
        #self.cluster_monospecies_count = 0  # private clusters (mono-species)
        #self.cluster_multispecies_ids = set()
        #self.cluster_multispecies_count = 0  # multispecies clusters
        #self.protein_monospecies_count = 0  # proteins in private clusters (mono-species)
        #self.protein_multispecies_count = 0  # proteins in private clusters (mono-species)
        #self.protein_singleton_count = 0   # singleton proteins (unclustered)

### more general RankInstanceObj()
'''
RankInstanceObj:
    - ProteomeObj, GroupObj, SpeciesObj, ...
    - rank id
    - members
    - clusters to which the at least one member of rank belongs
        - by type (singleton, monoton, multiton) : list of clusters
        - by percentage (only monotons, multitons) : list of clusters
        - by counts :
    - proteins belonging to rank in clusters of type
        - by type (singleton, monoton, multiton) : list of proteins
    -

'''
class RankInstanceObj():
    def __init__(self, rank_id):
        self.rank_id = rank_id
        self.members = members

class ProteinObj():
    def __init__(self, proteinID, length, clusterID):
        self.proteinID = proteinID
        self.length = length
        self.clusterID = clusterID

class ClusterObj():
    def __init__(self, clusterID, proteins):
        self.clusterID = clusterID
        self.proteins = proteins
        self.type_by_rank = {}

class CategoryObj():
    def __init__(self, members_by_rank, ranks_by_proteome):
        self.ranks = set(members_by_rank.keys())                                        # rank_keys
        self.members_by_rank = members_by_rank                                          # rank_key : rank_values
        self.member_count_by_rank = {x : len(v) for x, v in members_by_rank.items()}    # rank_key : members count
        self.ranks_by_proteome = ranks_by_proteome                                      # proteome : rank_key : rank_value

class ClusterCountScheme():
    def __init__(self, categoryObj):
        self.members_at_rank = {}
        self.member_count_at_rank = {}
        self.member_perc_at_rank = {}
        self.cluster_type_at_rank = {}

def parse_classification(species_classification_f):
    ranks_by_proteome = {}
    members_by_rank = {}
    header = []
    with open(species_classification_f) as fh:
        for l in fh:
            if l.startswith("#"):
                header = [x.strip() for x in l.lstrip("#").rstrip("\n").split()]
                members_by_rank = {rank : set() for rank in header}
            else:
                line = l.rstrip("\n").split()
                for idx, item in enumerate(line):
                    members_by_rank[header[idx]].add(item)
                    if not line[0] in ranks_by_proteome:
                        ranks_by_proteome[line[0]] = {}
                    ranks_by_proteome[line[0]][header[idx]] = item
    return members_by_rank, ranks_by_proteome


if __name__ == "__main__":

    species_ids_f, species_classification_f, groups_fs = '','',''
    try:
        species_ids_f = sys.argv[1]
        species_classification_f = sys.argv[2]
        groups_fs = sys.argv[3:] # one or more, format "OrthologousGroups_INFLATIONVALUE.txt"
    except:
        sys.exit("./kinfin.py SPECIESID_FILE SPECIESCLASSIFICATION_FILE GROUPS_FILE")

    PERCENTAGES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.55, 0.60, 0.65, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0]

    # Get all info from species_classification_f
    members_by_rank, ranks_by_proteome = parse_classification(species_classification_f)
    # Initialise CategoryObj with the info
    categoryObj = CategoryObj(members_by_rank, ranks_by_proteome)

    # Initialise DataObj
    dataObj = DataObj()
    # Add proteome-info to DataObj
    dataObj.add_proteomeObjs(species_ids_f)
    # Add category-info to DataObj (CategoryObj.ranks_by_proteome())
    dataObj.add_ranks_to_proteomeObjs(categoryObj)


    dataObj.setup_dirs()
    dataObj.add_orthologous_groups(groups_fs)
    #dataObj.update_clusterObjs()

    #for proteomeObj in dataObj.yield_proteomeObjs():
    #    print proteomeObj.__dict__
    #for inflation_value, clusterObj in dataObj.yield_clusterObj():
    #    print inflation_value, clusterObj.__dict__

    sys.exit("done")
    # for filtering clusters that have a mean min length of MIN_MEAN_LENGTH_OF_PROTEIN
    # - provide fasta length file
    MIN_MEAN_LENGTH_OF_PROTEIN = 40 # and less than one ASSEMBLY



    results = generate_dataObjs(species_ids_f, groups_fs)

    output(results)
