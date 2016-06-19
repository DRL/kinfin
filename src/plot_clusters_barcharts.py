#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import sys
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib as mat
mat.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("talk")
#sns.set(font='serif')
sns.set(style="whitegrid")

def plot_counts(infile, format):

    data = pd.read_csv(infile, sep="\t")

    f, ax = plt.subplots(figsize=(24, 12))
    sns.set_color_codes("pastel")
    cluster_singleton_count = data["cluster_singleton_count"] + data["cluster_monospecies_count"] + data["cluster_multispecies_count"]
    sns.barplot(x="inflation_value", y=cluster_singleton_count, data=data, label="Singleton clusters", color="g")
    cluster_monospecies_count = data["cluster_multispecies_count"] + data["cluster_monospecies_count"]
    sns.barplot(x="inflation_value", y=cluster_monospecies_count, data=data, label="Mono-species clusters", color="r")
    sns.barplot(x="inflation_value", y="cluster_multispecies_count", data=data, label="Multi-species clusters", color="b")
    ax.legend(ncol=3, loc="upper center", frameon=True)
    f.tight_layout()
    f.savefig("barplot.clusters.count." + format)

    g, ax = plt.subplots(figsize=(24, 12))
    protein_singleton_count = data["protein_singleton_count"] + data["protein_monospecies_count"] + data["protein_multispecies_count"]
    sns.barplot(x="inflation_value", y=protein_singleton_count, data=data, label="Singleton proteins", color="g")
    protein_monospecies_count = data["protein_multispecies_count"] + data["protein_monospecies_count"]
    sns.barplot(x="inflation_value", y=protein_monospecies_count, data=data, label="Mono-species proteins", color="r")
    sns.barplot(x="inflation_value", y="protein_multispecies_count", data=data, label="Multi-species proteins", color="b")
    ax.legend(ncol=3, loc="upper center", frameon=True)
    g.tight_layout()
    g.savefig("barplot.clusters.proteins.count." + format)




def plot_fraction(infile, format):

    data = pd.read_csv(infile, sep="\t")
    f, ax = plt.subplots(figsize=(24, 12))
    sns.set_color_codes("pastel")
    cluster_singleton_fraction = data["cluster_singleton_count"]/data["cluster_count"] + data["cluster_monospecies_count"]/data["cluster_count"] + data["cluster_multispecies_count"]/data["cluster_count"]
    sns.barplot(x="inflation_value", y=cluster_singleton_fraction, data=data, label="Singleton clusters", color="g")
    cluster_monospecies_fraction = data["cluster_multispecies_count"]/data["cluster_count"] + data["cluster_monospecies_count"]/data["cluster_count"]
    sns.barplot(x="inflation_value", y=cluster_monospecies_fraction, data=data, label="Mono-species clusters", color="r")
    cluster_multispecies_fraction = data["cluster_multispecies_count"]/data["cluster_count"]
    sns.barplot(x="inflation_value", y=cluster_multispecies_fraction, data=data, label="Multi-species clusters", color="b")
    ax.legend(ncol=3, loc="upper center", frameon=True)
    f.tight_layout()
    f.savefig("barplot.clusters.fraction." + format)

    g, ax = plt.subplots(figsize=(24, 12))
    protein_singleton_fraction = data["protein_singleton_count"]/data["protein_count"] + data["protein_monospecies_count"]/data["protein_count"] + data["protein_multispecies_count"]/data["protein_count"]
    sns.barplot(x="inflation_value", y=protein_singleton_fraction, data=data, label="Singleton proteins", color="g")
    protein_monospecies_fraction = data["protein_multispecies_count"]/data["protein_count"] + data["protein_monospecies_count"]/data["protein_count"]
    sns.barplot(x="inflation_value", y=protein_monospecies_fraction, data=data, label="Mono-species proteins", color="r")
    protein_multispecies_fraction = data["protein_multispecies_count"]/data["protein_count"]
    sns.barplot(x="inflation_value", y=protein_multispecies_fraction, data=data, label="Multi-species proteins", color="b")
    ax.legend(ncol=3, loc="upper center", frameon=True)
    g.tight_layout()
    g.savefig("barplot.clusters.proteins.fraction." + format)

if __name__ == "__main__":

    try:
        infile = sys.argv[1]
        format = sys.argv[2]
    except:
        sys.exit("./plot_groups_barchart.py INFILE FORMAT")

    plot_counts(infile, format)
    plot_fraction(infile, format)
