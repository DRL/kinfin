#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import sys
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.colors as colors
import matplotlib as mat
mat.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("talk")
#sns.set(font='serif')
sns.set(style="whitegrid")


def plot(out_f, files):
    f, ax = plt.subplots(figsize=(24, 12))
    sns.set_color_codes("pastel")
    for idx, infile in enumerate(files):
        colour = plt.cm.Paired(idx/len(files))
        print infile, colors.rgb2hex(colour)
        x = 0
        x_values = []
        y_values = []
        with open(infile) as fh:
            for line in fh:
                y_values.append(line.rstrip("\n"))
                x_values.append(x)
                x +=1
        ax.plot(x_values, y_values, '-o', label=infile)
        ax.set_xlim([0, len(x_values)])
        ax.set_ylim([0, 1.1])
        f.tight_layout()
    out_f = "%s.coverage.curve.png" % out_f
    ax.legend(ncol=3, loc="upper center", frameon=True)
    f.tight_layout()
    f.savefig(out_f, format="png")
        #ax = plt.errorbar(data["sample_size"], data["clusters_data_all_mean"], yerr=[data["clusters_data_all_max"], data["clusters_data_all_min"]])


    #cluster_monospecies_count = data["cluster_multispecies_count"] + data["cluster_monospecies_count"]
    #sns.barplot(x="inflation_value", y=cluster_monospecies_count, data=data, label="Mono-species clusters", color="r")
    #sns.barplot(x="inflation_value", y="cluster_multispecies_count", data=data, label="Multi-species clusters", color="b")
    #ax.legend(ncol=3, loc="upper center", frameon=True)


#
    #g, ax = plt.subplots(figsize=(24, 12))
#
    #protein_singleton_count = data["protein_singleton_count"] + data["protein_monospecies_count"] + data["protein_multispecies_count"]
    #sns.barplot(x="inflation_value", y=protein_singleton_count, data=data, label="Singleton proteins", color="g")
    #protein_monospecies_count = data["protein_multispecies_count"] + data["protein_monospecies_count"]
    #sns.barplot(x="inflation_value", y=protein_monospecies_count, data=data, label="Mono-species proteins", color="r")
    #sns.barplot(x="inflation_value", y="protein_multispecies_count", data=data, label="Multi-species proteins", color="b")
    ##sns.despine(left=False, offset=10, trim=True)
    #ax.legend(ncol=3, loc="upper center", frameon=True)
    #ax3.set_ylabel("Count")
    #ax3.set_xlabel("")
#

    #g.tight_layout()
    #g.savefig("barplot.clusters.proteins.svg")




    #perc_count_monospecies = data["count_monospecies"]/data["count_clusters"]
    #clusters = data["count_clusters"]
    #perc_count_singletons = data["count_singletons"]/data["count_clusters"]
    #inflation_values = data["#inflation_value"]
    #dataframe = pd.concat([clusters, perc_count_monospecies, perc_count_singletons, inflation_values], axis=1, keys=['clusters', 'mono-species clusters (%)', 'single-protein clusters (%)', 'Inflation value'])
    #f = sns.lmplot('mono-species clusters (%)', 'clusters', data=dataframe, hue='Inflation value', fit_reg=False)
    #f.add_legend()
    #f.savefig("mono_vs_clusters.png")
    #g = sns.lmplot('mono-species clusters (%)', 'single-protein clusters (%)', data=dataframe, hue='Inflation value', fit_reg=False)
    #g.add_legend()
    #g.savefig("mono_vs_single.png")
    #h = sns.lmplot('single-protein clusters (%)', 'clusters',  data=dataframe, hue='Inflation value', fit_reg=False)
    #h.add_legend()
    #h.savefig("single_vs_clusters.png")

if __name__ == "__main__":

    try:
        out_f = sys.argv[1]
        files = sys.argv[2:]
    except:
        sys.exit("./seaborn_histogram.py INFILE")

    plot(out_f, files)
