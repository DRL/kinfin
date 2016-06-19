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


def plot(files):
    f, ax = plt.subplots(figsize=(24, 12))
    sns.set_color_codes("pastel")
    labels = ["1.5 => 2.0", "2.0 => 2.5", "2.5 => 3.0", "3.0 => 3.5", "3.5 => 4.0", "4.0 => 4.5", "4.5 => 5.0"]
    for idx, infile in enumerate(files):
        colour = plt.cm.Paired(idx/len(files))
        print infile, colors.rgb2hex(colour)
        data = pd.read_csv(infile, sep="\t")
        #print data.head()


        #ax = plt.errorbar(data["sample_size"], data["clusters_data_all_mean"], yerr=[data["clusters_data_all_max"], data["clusters_data_all_min"]])

        ax = sns.distplot(data, kde=False, bins=100, label=labels[idx], hist_kws={"histtype": "barstacked", "linewidth": 1, "color": colour})
    #cluster_monospecies_count = data["cluster_multispecies_count"] + data["cluster_monospecies_count"]
    #sns.barplot(x="inflation_value", y=cluster_monospecies_count, data=data, label="Mono-species clusters", color="r")
    #sns.barplot(x="inflation_value", y="cluster_multispecies_count", data=data, label="Multi-species clusters", color="b")
    ax.legend(ncol=1, loc="upper center", frameon=True)
    f.get_axes()[0].set_yscale('log')
    f.tight_layout()
    f.savefig("jaccard.histogram.png")
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
        files = sys.argv[1:]
    except:
        sys.exit("./seaborn_histogram.py INFILES")

    plot(files)
