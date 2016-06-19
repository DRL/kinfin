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
from matplotlib.colors import colorConverter
import seaborn as sns
sns.set_context("talk")
#sns.set(font='serif')
sns.set(style="whitegrid")

def plot_counts(infile, format):

    data = pd.read_csv(infile, sep="\t")
    data = data.sort_values(['species_class', 'species_id'], ascending=[True, True])
    f, ax = plt.subplots(figsize=(24, 12))
    species_protein_singleton_fraction = data["species_protein_singleton_fraction"] * data["species_protein_count"] + data["species_protein_monospecies_fraction"] * data["species_protein_count"] + data["species_protein_multispecies_fraction"] * data["species_protein_count"]
    sns.barplot(x="species_id", y=species_protein_singleton_fraction, data=data,
            label="singletons", color=colorConverter.to_rgb('#5599FF'))
    species_protein_monospecies_fraction = data["species_protein_monospecies_fraction"] * data["species_protein_count"] + data["species_protein_multispecies_fraction"] * data["species_protein_count"]
    sns.barplot(x="species_id", y=species_protein_monospecies_fraction,  data=data,
            label="monospecies", color=colorConverter.to_rgb('#FF9955'))
    species_protein_multispecies_fraction = data["species_protein_multispecies_fraction"] * data["species_protein_count"]
    sns.barplot(x="species_id", y=species_protein_multispecies_fraction,  data=data,
            label="multispecies", color=colorConverter.to_rgb('#FFCC00'))
    ax.legend(ncol=1, loc="upper right", frameon=True)
    ax.set_ylabel("Count of proteins in clusters")
    ax.set_xlabel("Inflation value")
    sns.despine(left=False, offset=10, trim=True)
    plt.xticks(rotation=90)
    f.tight_layout()

    f.savefig("barplot.species.counts_proteins." + format)

def plot_fraction(infile, format):
    data = pd.read_csv(infile, sep="\t")
    data = data.sort_values(['species_class', 'species_id'], ascending=[True, True])
    f, ax = plt.subplots(figsize=(24, 12))
    sns.set_color_codes("pastel")

    species_protein_singleton_fraction = data["species_protein_singleton_fraction"] + data["species_protein_monospecies_fraction"] + data["species_protein_multispecies_fraction"]
    species_protein_singleton_fraction = np.ones(121) # hack until kinfin output fixed
    sns.barplot(x="species_id", y=species_protein_singleton_fraction, data=data,
            label="singletons", color=colorConverter.to_rgb('#5599FF'))
    species_protein_monospecies_fraction = data["species_protein_monospecies_fraction"] + data["species_protein_multispecies_fraction"]
    sns.barplot(x="species_id", y=species_protein_monospecies_fraction,  data=data,
            label="monospecies", color=colorConverter.to_rgb('#FF9955'))
    species_protein_multispecies_fraction = data["species_protein_multispecies_fraction"]
    sns.barplot(x="species_id", y=species_protein_multispecies_fraction,  data=data,
            label="multispecies", color=colorConverter.to_rgb('#FFCC00'))
    ax.legend(ncol=1, loc="upper right", frameon=True)
    ax.set_ylabel("Percentage of proteins in clusters")
    ax.set_xlabel("Inflation value")
    sns.despine(left=False, offset=10, trim=True)
    plt.xticks(rotation=90)
    f.tight_layout()

    f.savefig("barplot.species.fraction_proteins." + format)

if __name__ == "__main__":

    try:
        infile = sys.argv[1]
        format = sys.argv[2]
    except:
        sys.exit("./plot.py INFILE FORMAT")

    plot_counts(infile, format)
    plot_fraction(infile, format)
