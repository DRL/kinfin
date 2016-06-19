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

def plot(infile1_5, infile2_0, infile2_5, infile3_0, infile3_5, infile4_0, infile4_5, infile5_0, format):

    #head -1 cluster_stats.txt > cluster_stats.1_5.txt ; awk '$1 == "I1.5"' cluster_stats.txt >> cluster_stats.1_5.txt; \
    #head -1 cluster_stats.txt > cluster_stats.2_0.txt ; awk '$1 == "I2.0"' cluster_stats.txt >> cluster_stats.2_0.txt; \
    #head -1 cluster_stats.txt > cluster_stats.2_5.txt ; awk '$1 == "I2.5"' cluster_stats.txt >> cluster_stats.2_5.txt; \
    #head -1 cluster_stats.txt > cluster_stats.3_0.txt ; awk '$1 == "I3.0"' cluster_stats.txt >> cluster_stats.3_0.txt; \
    #head -1 cluster_stats.txt > cluster_stats.3_5.txt ; awk '$1 == "I3.5"' cluster_stats.txt >> cluster_stats.3_5.txt; \
    #head -1 cluster_stats.txt > cluster_stats.4_0.txt ; awk '$1 == "I4.0"' cluster_stats.txt >> cluster_stats.4_0.txt; \
    #head -1 cluster_stats.txt > cluster_stats.4_5.txt ; awk '$1 == "I4.5"' cluster_stats.txt >> cluster_stats.4_5.txt; \
    #head -1 cluster_stats.txt > cluster_stats.5_0.txt ; awk '$1 == "I5.0"' cluster_stats.txt >> cluster_stats.5_0.txt

    data1_5 = pd.read_csv(infile1_5, sep="\t")
    data2_0 = pd.read_csv(infile2_0, sep="\t")
    data2_5 = pd.read_csv(infile2_5, sep="\t")
    data3_0 = pd.read_csv(infile3_0, sep="\t")
    data3_5 = pd.read_csv(infile3_5, sep="\t")
    data4_0 = pd.read_csv(infile4_0, sep="\t")
    data4_5 = pd.read_csv(infile4_5, sep="\t")
    data5_0 = pd.read_csv(infile5_0, sep="\t")

    f, ax = plt.subplots(figsize=(24, 12))
    sns.set_color_codes("pastel")

    ax = sns.distplot(data1_5["cluster_protein_count"], kde=False, bins=100, label="I1_5", hist_kws={"histtype": "step", "linewidth": 2})
    ax = sns.distplot(data2_0["cluster_protein_count"], kde=False, bins=100, label="I2_0", hist_kws={"histtype": "step", "linewidth": 2})
    ax = sns.distplot(data2_5["cluster_protein_count"], kde=False, bins=100, label="I2_5", hist_kws={"histtype": "step", "linewidth": 2})
    ax = sns.distplot(data3_0["cluster_protein_count"], kde=False, bins=100, label="I3_0", hist_kws={"histtype": "step", "linewidth": 2})
    ax = sns.distplot(data3_5["cluster_protein_count"], kde=False, bins=100, label="I3_5", hist_kws={"histtype": "step", "linewidth": 2})
    ax = sns.distplot(data4_0["cluster_protein_count"], kde=False, bins=100, label="I4_0", hist_kws={"histtype": "step", "linewidth": 2})
    ax = sns.distplot(data4_5["cluster_protein_count"], kde=False, bins=100, label="I4_5", hist_kws={"histtype": "step", "linewidth": 2})
    ax = sns.distplot(data5_0["cluster_protein_count"], kde=False, bins=100, label="I5_0", hist_kws={"histtype": "step", "linewidth": 2})

    ax.legend(ncol=1, loc="upper right", frameon=True)
    ax.set_ylabel("Count")
    ax.set_xlabel("")
    ax.set_yscale('log')
    ax.set_ylim(0.1,)
    sns.despine(left=False, offset=10, trim=True)
    plt.xticks(rotation=90)
    f.tight_layout()

    f.savefig("histogram_of_protein_count_by_inflation_value.png")
    f.savefig("histogram_of_protein_count_by_inflation_value.svg")

if __name__ == "__main__":
    try:
        infile1_5, infile2_0, infile2_5, infile3_0, infile3_5, infile4_0, infile4_5, infile5_0 = sys.argv[1:9]
        format = sys.argv[2]
    except:
        sys.exit("./seaborn_histogram.py INFILE")
    print infile1_5
    print infile2_0
    print infile2_5
    print infile3_0
    print infile3_5
    print infile4_0
    print infile4_5
    print infile5_0
    plot(infile1_5, infile2_0, infile2_5, infile3_0, infile3_5, infile4_0, infile4_5, infile5_0, format)
