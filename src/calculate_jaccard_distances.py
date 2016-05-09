#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import sys
import os
import numpy as np
import itertools

class CountObj:
    def __init__(self, inflation_value_combination):
        self.inflation_value_combination = inflation_value_combination
        self.count_of_clusters_A = 0
        self.count_of_clusters_B = 0
        self.count_of_clusters_identical = 0
        self.count_of_proteins_identical = 0
        self.count_of_clusters_different_hierarchical = 0
        self.count_of_proteins_different_hierarchical = 0
        self.count_of_clusters_different_not_hierarchical = 0
        self.count_of_proteins_different_not_hierarchical = 0
        self.count_of_singletons_identical = 0
        self.count_of_singletons_different = 0

def calculate_jaccards(clusters, cluster_ids, inflation_value_A, inflation_value_B):
    print "calculating"
    jaccard_distances_of_combinations = {}
    #inflation_values = [("I1.5", "I2.0"), ("I2.0", "I2.5"), ("I2.5", "I3.5"), ("I3.5", "I4.0"), ("I4.0", "I5.0")]
    inflation_values = [(inflation_value_A, inflation_value_B)]
    counts = {}
    for inflation_value_combination in inflation_values:

        print
        print inflation_value_combination

        list_of_jaccard_distances = []
        A, B = inflation_value_combination

        countObj = CountObj(inflation_value_combination)
        countObj.count_of_clusters_A = len(clusters[A])
        countObj.count_of_clusters_B = len(clusters[B])

        clusters_in_A_different = []

        for idx, cluster_A in enumerate(clusters[A]):
            # Progress
            #if (idx % 1000) == 0:
            #    sys.stdout.write('\r')
            #    print '{:.2%}'.format(idx/countObj.count_of_clusters_A),
            #    sys.stdout.flush()
            # setup
            cardinality_cluster_A = len(cluster_A)
            sum_of_jaccard_distances = 0
            # cluster is multiton
            clusters_identical = 0
            if cardinality_cluster_A > 1:
                for cluster_B in clusters[B]:
                    if sum_of_jaccard_distances == 1.0:
                        break
                    intersection = cluster_A & cluster_B
                    if (intersection):
                        cardinality_cluster_B = len(cluster_B)
                        cardinality_intersection = len(intersection)
                        # clusters are identical
                        if cardinality_intersection == cardinality_cluster_A:
                            countObj.count_of_clusters_identical += 1
                            countObj.count_of_proteins_identical += cardinality_intersection
                            sum_of_jaccard_distances += 1.0
                            clusters_identical = 1
                            break
                        # clusters are different but intersect in more than one protein
                        else:
                            if cardinality_intersection == cardinality_cluster_B:
                                #countObj.count_of_clusters_different_hierarchical += 1
                                countObj.count_of_proteins_different_hierarchical += cardinality_intersection
                            else:
                                #countObj.count_of_clusters_different_not_hierarchical += 1
                                countObj.count_of_proteins_different_not_hierarchical += cardinality_intersection
                            sum_of_jaccard_distances += cardinality_intersection/(cardinality_cluster_A + cardinality_cluster_B - cardinality_intersection)
                if sum_of_jaccard_distances == 1.0:
                    if not clusters_identical:
                        countObj.count_of_clusters_different_hierarchical += 1
                else:
                    countObj.count_of_clusters_different_not_hierarchical += 1
                    clusters_in_A_different.append(cluster_ids[inflation_value_combination[0]][idx])
            # cluster is singleton
            else:
                for cluster_B in clusters[B]:
                    intersection = cluster_A & cluster_B
                    if (intersection):
                        cardinality_cluster_B = len(cluster_B)
                        if cardinality_cluster_B == 1:
                            sum_of_jaccard_distances += 1.0
                            countObj.count_of_clusters_identical += 1
                            countObj.count_of_proteins_identical += 1
                            countObj.count_of_singletons_identical += 1
                            break
                        elif cardinality_cluster_B > 1:
                            sum_of_jaccard_distances += 1/(1 + cardinality_cluster_B - 1)
                            countObj.count_of_clusters_different_not_hierarchical += 1
                            countObj.count_of_proteins_different_not_hierarchical += 1
                            countObj.count_of_singletons_different += 1
                            clusters_in_A_different.append(cluster_ids[inflation_value_combination[0]][idx])
                            break
                        else:
                            sys.exit("ERROR3")
            list_of_jaccard_distances.append((cardinality_cluster_A, sum_of_jaccard_distances))
        counts[inflation_value_combination] = countObj
        with open("jaccard_" + "to".join(inflation_value_combination) + ".txt", "w") as fh:
            for jaccard_distance in list_of_jaccard_distances:
                cardinality_cluster_A = jaccard_distance[0]
                sum_of_jaccard_distance = jaccard_distance[1]
                fh.write("%s\t%s\n" % (cardinality_cluster_A, sum_of_jaccard_distance))
        with open("jaccard_not_hierarchical_" + "to".join(inflation_value_combination) + ".txt", "w") as fh:
            for cluster_id in clusters_in_A_different:
                fh.write("%s\n" % cluster_id)

        with open("jaccard_stats." + "to".join(inflation_value_combination) + ".txt", "w") as fh:
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (\
                "inflation_value_combination", \
                "count_of_clusters_A", \
                "count_of_clusters_B", \
                "count_of_clusters_identical", \
                "count_of_proteins_identical", \
                "count_of_clusters_different_hierarchical", \
                "count_of_proteins_different_hierarchical", \
                "count_of_clusters_different_not_hierarchical", \
                "count_of_proteins_different_not_hierarchical", \
                "count_of_singletons_identical", \
                "count_of_singletons_different" \
                ))
            for inflation_value_combination in sorted(counts):
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (\
                    "to".join(inflation_value_combination), \
                    counts[inflation_value_combination].count_of_clusters_A, \
                    counts[inflation_value_combination].count_of_clusters_B, \
                    counts[inflation_value_combination].count_of_clusters_identical, \
                    counts[inflation_value_combination].count_of_proteins_identical, \
                    counts[inflation_value_combination].count_of_clusters_different_hierarchical, \
                    counts[inflation_value_combination].count_of_proteins_different_hierarchical, \
                    counts[inflation_value_combination].count_of_clusters_different_not_hierarchical, \
                    counts[inflation_value_combination].count_of_proteins_different_not_hierarchical, \
                    counts[inflation_value_combination].count_of_singletons_identical, \
                    counts[inflation_value_combination].count_of_singletons_different \
                    ))

def parse_cluster_stats(cluster_stats_f, inflation_value_A, inflation_value_B):
    clusters = {}
    cluster_ids = {}
    print "parsing"
    with open(cluster_stats_f) as fh:
        for l in fh:
            if not l.startswith("#"):
                line = l.rstrip("\n").split()
                inflation_value = line[0]
                proteins = set(line[9].split(","))
                cluster_id = line[1]
                #proteins = set(line[7].split(","))
                if inflation_value == inflation_value_A or inflation_value == inflation_value_B:
                    if not inflation_value in clusters:
                        print inflation_value
                        clusters[inflation_value] = []
                        cluster_ids[inflation_value] = []
                    clusters[inflation_value].append(proteins)
                    cluster_ids[inflation_value].append(cluster_id)

    return clusters, cluster_ids

#def output(jaccard_distances_of_combinations):
#    for combination, jaccard_distances in jaccard_distances_of_combinations.items():
#        with open(combination + ".txt", "w") as fh:
#            for jaccard_distance in jaccard_distances:
#                fh.write("%s\n" % jaccard_distance)


if __name__ == "__main__":

    try:
        cluster_stats_f = sys.argv[1]
        inflation_value_A = sys.argv[2]
        inflation_value_B = sys.argv[3]
    except:
        sys.exit("./calculate_jaccards.py cluster_stats.txt inflation_value_A inflation_value_B")

    clusters, cluster_ids = parse_cluster_stats(cluster_stats_f, inflation_value_A, inflation_value_B)
    calculate_jaccards(clusters, cluster_ids, inflation_value_A, inflation_value_B)
    #output(jaccard_distances_of_combinations)
