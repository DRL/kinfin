#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: get_protein_cluster_count.py     -c <FILE> -e <STRING>
                                        [-o <OUTPREFIX>]
                                        [-h|--help]

    Options:
        -h --help                           show this

        -c <FILE>                           Counts of clusters file
        -e, --exclude <STRING>              Exclude taxon [default: "CELEG"]
        -o, --outprefix <OUTPREFIX>         Prefix for output files
"""

from __future__ import division
from docopt import docopt
import sys
from os.path import basename, isfile, abspath, splitext, join, exists
import numpy as np
from scipy import arange, stats
import scipy.cluster.hierarchy as sch


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
import pylab
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
FONTSIZE = 36
plt.rcParams['font.size'] = FONTSIZE
plt.rcParams['axes.labelsize'] = FONTSIZE
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = FONTSIZE
plt.rcParams['xtick.labelsize'] = FONTSIZE-2
plt.rcParams['ytick.labelsize'] = FONTSIZE-2
plt.rcParams['legend.fontsize'] = FONTSIZE
plt.rcParams['figure.titlesize'] = FONTSIZE+2

from math import *
from random import *

# function x=randht(n, varargin)

# RANDHT generates n observations distributed as some continous heavy-
# tailed distribution. Options are power law, log-normal, stretched
# exponential, power law with cutoff, and exponential. Can specify lower
# cutoff, if desired.
#
#    Example:
#       x = randht(10000,'powerlaw',alpha);
#       x = randht(10000,'xmin',xmin,'powerlaw',alpha);
#       x = randht(10000,'cutoff',alpha, lambda);
#       x = randht(10000,'exponential',lambda);
#       x = randht(10000,'lognormal',mu,sigma);
#       x = randht(10000,'stretched',lambda,beta);
#
#    See also PLFIT, PLVAR, PLPVA
#
#    Source: http://www.santafe.edu/~aaronc/powerlaws/


# Version 1.0.2 (2008 April)
# Copyright (C) 2007 Aaron Clauset (Santa Fe Institute)

# Ported to python by Joel Ornstein (2011 August)
# (joel_ornstein@hmc.edu)

# Distributed under GPL 2.0
# http://www.gnu.org/copyleft/gpl.html
# RANDHT comes with ABSOLUTELY NO WARRANTY
#
# Notes:
#
def randht(n,*varargin):
    Type   = '';
    xmin   = 1;
    alpha  = 2.5;
    beta   = 1;
    Lambda = 1;
    mu     = 1;
    sigma  = 1;


    # parse command-line parameters; trap for bad input
    i=0;
    while i<len(varargin):
        argok = 1;
        if type(varargin[i])==str:
            if varargin[i] == 'xmin':
                xmin = varargin[i+1]
                i = i + 1
            elif varargin[i] == 'powerlaw':
                Type = 'PL'
                alpha  = varargin[i+1]
                i = i + 1
            elif varargin[i] == 'cutoff':
                Type = 'PC';
                alpha  = varargin[i+1]
                Lambda = varargin[i+2]
                i = i + 2
            elif varargin[i] == 'exponential':
                Type = 'EX'
                Lambda = varargin[i+1]
                i = i + 1
            elif varargin[i] == 'lognormal':
                Type = 'LN';
                mu = varargin[i+1]
                sigma = varargin[i+2]
                i = i + 2
            elif varargin[i] == 'stretched':
                Type = 'ST'
                Lambda = varargin[i+1]
                beta = varargin[i+2]
                i = i + 2
            else: argok=0


        if not argok:
            print '(RANDHT) Ignoring invalid argument #' ,i+1

        i = i+1

    if n<1:
        print '(RANDHT) Error: invalid ''n'' argument; using default.\n'
        n = 10000;

    if xmin < 1:
        print '(RANDHT) Error: invalid ''xmin'' argument; using default.\n'
        xmin = 1;




    x=[]
    if Type == 'EX':
        x=[]
        for i in range(n):
            x.append(xmin - (1./Lambda)*log(1-random()))
    elif Type == 'LN':
        y=[]
        for i in range(10*n):
            y.append(exp(mu+sigma*normalvariate(0,1)))

        while True:
            y= filter(lambda X:X>=xmin,y)
            q = len(y)-n;
            if q==0: break

            if q>0:
                r = range(len(y));
                shuffle(r)
                ytemp = []
                for j in range(len(y)):
                    if j not in r[0:q]:
                        ytemp.append(y[j])
                y=ytemp
                break
            if (q<0):
                for j in range(10*n):
                    y.append(exp(mu+sigma*normalvariate(0,1)))

        x = y

    elif Type =='ST':
        x=[]
        for i in range(n):
            x.append(pow(pow(xmin,beta) - (1./Lambda)*log(1.-random()),(1./beta)))
    elif Type == 'PC':

        x = []
        y=[]
        for i in range(10*n):
            y.append(xmin - (1./Lambda)*log(1.-random()))
        while True:
            ytemp=[]
            for i in range(10*n):
                if random()<pow(y[i]/float(xmin),-alpha):ytemp.append(y[i])
            y = ytemp
            x = x+y
            q = len(x)-n
            if q==0: break;

            if (q>0):
                r = range(len(x))
                shuffle(r)

                xtemp = []
                for j in range(len(x)):
                    if j not in r[0:q]:
                        xtemp.append(x[j])
                x=xtemp
                break;

            if (q<0):
                y=[]
                for j in range(10*n):
                    y.append(xmin - (1./Lambda)*log(1.-random()))


    else:
        x=[]
        for i in range(n):
            x.append(xmin*pow(1.-random(),-1./(alpha-1.)))

    return x

def parse_counts(counts_of_clusters_f):
    counts_ALL = {}
    counts_ONLY_INCLUDE = {}
    EXCLUDE_idxs = []
    INCLUDE_idxs = []
    with open(counts_of_clusters_f) as in_fh:
        for line in in_fh:
            temp = line.rstrip("\n").split()
            if line.startswith("#"):
                for idx, element in enumerate(temp):
                    if element in EXCLUDE:
                        print "[+] - Excluding %s" % element
                        EXCLUDE_idxs.append(idx)
                if not EXCLUDE_idxs:
                    sys.exit("[ERROR] - %s could not be found in header" % EXCLUDE)
                INCLUDE_idxs = [x for x in range(1, len(temp)) if x not in EXCLUDE_idxs]
                for INCLUDE_idx in INCLUDE_idxs:
                    print "[+] - Including %s" % temp[INCLUDE_idx]
            else:
                counts_cluster = int(temp[0])
                counts_protein = sum([int(x) for x in temp[1:]])
                if sum([int(temp[EXCLUDE_idx]) for EXCLUDE_idx in EXCLUDE_idxs]) == 0:
                    counts_ONLY_INCLUDE[counts_protein] = counts_ONLY_INCLUDE.get(counts_protein, 0) + counts_cluster
                counts_ALL[counts_protein] = counts_ALL.get(counts_protein, 0) + counts_cluster

#    if (prefix):
#        counts_with_f = "%s.%s" % (prefix, "counts_with.txt")
#        counts_without_f = "%s.%s" % (prefix, "counts_without.txt")
#    else:
#        counts_with_f = "%s" % ("counts_with.txt")
#        counts_without_f = "%s" % ("counts_without.txt")
#    with open(counts_with_f, "w") as counts_with_fh:
#        for count_WITH_EXCLUDE in sorted(counts_WITH_EXCLUDE):
#            counts_with_fh.write("%s\t%s\n" % (count_WITH_EXCLUDE, counts_WITH_EXCLUDE[count_WITH_EXCLUDE]))
#    with open(counts_without_f, "w") as counts_without_fh:
#        for count_WITHOUT_EXCLUDE in sorted(counts_WITHOUT_EXCLUDE):
#            counts_without_fh.write("%s\t%s\n" % (count_WITHOUT_EXCLUDE, counts_WITHOUT_EXCLUDE[count_WITHOUT_EXCLUDE]))
    return counts_ALL, counts_ONLY_INCLUDE

def plot(counts_ALL, counts_ONLY_INCLUDE):
    f, ax = plt.subplots(figsize=(24,18))
    counts_protein_values = []
    counts_family_values = []
    legend_handles = []
    legend_labels = []
    for count_ALL in sorted(counts_ALL):
        counts_family_values.append(counts_ALL[count_ALL])
        counts_protein_values.append(count_ALL)
    label = "present in any Metazoan proteome"
    plot = ax.plot(counts_protein_values, counts_family_values, 'o', markersize=20, label=label, markeredgecolor = 'none')
    #legend_handles.append(Line2D([0], [0], linestyle="none", marker="o", alpha=1, markersize=10))
    #legend_labels.append(label)

    counts_protein_values = []
    counts_family_values = []
    for count_ONLY_INCLUDE in sorted(counts_ONLY_INCLUDE):
        counts_family_values.append(counts_ONLY_INCLUDE[count_ONLY_INCLUDE])
        counts_protein_values.append(count_ONLY_INCLUDE)
    #label = "present in any Nematode species (absent in %s)" % ",\n".join(EXCLUDE)
    label = "only present in Nematode proteomes"
    plot = ax.plot(counts_protein_values, counts_family_values, 'd', markersize=20, label=label, markeredgecolor = 'none')

    #legend_handles.append(Line2D([0], [0], linestyle="none", marker="o", alpha=1, markersize=10))
    #legend_labels.append(label)
    ax.legend(ncol=1, numpoints=1, loc="upper right", frameon=True)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.margins(0.8)
    #minorLocator = MultipleLocator(100)
    #ax.xaxis.set_major_locator(majorLocator)
    #ax.yaxis.set_major_locator(majorLocator)
    #plt.gca().set_ylim(bottom=0.8, top=1000)
    plt.gca().set_ylim(bottom=0.8)
    #plt.gca().set_xlim(left=0.8, right=100)
    plt.gca().set_xlim(left=0.8)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_ylabel("Number of families")
    ax.set_xlabel("Size of family (Number of proteins)")
    f.tight_layout()
    ax.grid(True, which="minor", color="lightgrey")
    ax.grid(True, which="major", color="grey")
    f.savefig("%s.plot.png" % prefix, format="png")


if __name__ == "__main__":
    __version__ = 0.1
    counts_of_clusters_f = None
    EXCLUDE = None
    prefix = ''
    args = docopt(__doc__)

    counts_of_clusters_f = args['-c']
    EXCLUDE = args['--exclude']
    prefix = args['--outprefix']
    if "," in EXCLUDE:
        EXCLUDE = EXCLUDE.split(",")
    else:
        EXCLUDE = [EXCLUDE]
    counts_ALL, counts_ONLY_INCLUDE = parse_counts(counts_of_clusters_f)

    plot(counts_ALL, counts_ONLY_INCLUDE)
