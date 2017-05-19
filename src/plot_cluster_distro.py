#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: ips_to_table.py      -i <FILE> [-o <STR>]
                            [--domain_sources <STRING>]
                            [--only_domain_go]
                            [--only_domain_ipr]
                            [-h|--help]

    Options:
        -h --help                               show this
        -i, --interproscan_f <FILE>             Interproscan file
        -o, --out_prefix <STR>                  Outprefix (default: )
        --domain_sources <STRING>               Collect information for those domain_sources [default: SignalP_EUK,Pfam]
                                                    - GO-terms and IPR-IDs domain counts are inferred for these
"""

from __future__ import division
import sys
from collections import Counter

import_errors = []
try:
    from docopt import docopt
except ImportError:
    import_errors.append("[ERROR] : Module \'Docopt\' was not found. Please install \'Docopt\' using \'pip install docopt\'")
try:
    import scipy
except ImportError:
    import_errors.append("[ERROR] : Module \'SciPy\' was not found. Please install \'SciPy\' using \'pip install scipy\'")
try:
    import matplotlib as mat
except ImportError:
    import_errors.append("[ERROR] : Module \'Matplotlib\' was not found. Please install \'Matplotlib\' using \'pip install matplotlib\'")
if import_errors:
    sys.exit("\n".join(import_errors))

import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

mat.use('agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import seaborn as sns
sns.set_context("talk")
sns.set(style="whitegrid")
sns.set_color_codes("pastel")
mat.rc('ytick', labelsize=20)
mat.rc('xtick', labelsize=20)
axis_font = {'size':'20'}
mat.rcParams.update({'font.size': 22})

def parse_clusters(cluster_f):
    cluster_proteome_ids_by_cluster_id = {}
    with open(cluster_f) as cluster_fh:
        for line in cluster_fh:
            if not line.startswith("#"):
                temp = line.rstrip("\n").split("\t")
                cluster_id = temp[0]
                print temp
                proteome_ids = frozenset(temp[12].split(","))
                if temp[1] == 'present' and temp[2] == 'specific':
                    cluster_proteome_ids_by_cluster_id[cluster_id] = proteome_ids

    print "[STATUS] %s (%s) : %s specific-present clusters" % (hypothesis, ALO, len(cluster_proteome_ids_by_cluster_id))
    return cluster_proteome_ids_by_cluster_id

def calculate_results(cluster_proteome_ids_by_cluster_id, ALO):
    non_TA_counts_by_TA_counts = {}
    for cluster_id, cluster_proteome_ids in cluster_proteome_ids_by_cluster_id.items():
        TA_count = len(TA.intersection(cluster_proteome_ids))
        print TA_count
        non_TA_count = None
        if ALO == "NENP":
            non_TA_count = len(NENP.intersection(cluster_proteome_ids))
        elif ALO == "ARON":
            non_TA_count = len(ARON.intersection(cluster_proteome_ids))
        else:
            sys.exit("[ERROR] - %s" % ALO)
        if not TA_count in non_TA_counts_by_TA_counts:
            non_TA_counts_by_TA_counts[TA_count] = []
        non_TA_counts_by_TA_counts[TA_count].append(non_TA_count)
    return {TA_count : Counter(non_TA_counts) for TA_count, non_TA_counts in non_TA_counts_by_TA_counts.items()}

def plot_counts(non_TA_counts_by_TA_counts):
    non_TA_counts_by_TA_counts_plot_f = "%s.%s_counts_by_TA_counts.%s" % (hypothesis, ALO, PLOT_FORMAT)
    f, ax = plt.subplots(figsize=FIGSIZE)
    idx = 0
    for TA_count, non_TA_counter in non_TA_counts_by_TA_counts.items():
        clusters_count = 0
        #colour = plt.cm.Paired(idx/len(non_TA_counts_by_TA_counts))
        colour = plt.cm.Paired(idx/5)
        x_values = []
        y_values = []
        for value, count in non_TA_counter.items():
            if value > 0:
                clusters_count += count
            x_values.append(value)
            y_values.append(count)
        x_array = np.array(x_values)
        y_array = np.array(y_values)
        ax.plot(x_array, y_array, '-o', c=colour, alpha=0.8, label="%s (%s clusters w/ >= 1 non-TA)" % (TA_count, clusters_count))
        idx += 1
    ax.set_xlabel("Count of Non-Tardigrade members (%s)" % ALO, **axis_font)
    ax.set_ylabel("Count of clusters", **axis_font)
    ax.set_yscale('log')
    #ax.set_xscale('log')
    plt.margins(0.8)
    plt.gca().set_ylim(bottom=0.8)
    plt.gca().set_xlim(left=-0.2, right=max_non_TA_count+1)
    #plt.gca().set_xlim(left=0.8, right=100)
    #plt.gca().set_xlim(left=0.8)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    #f.tight_layout()
    ax.grid(True, linewidth=0.5, which="minor", color="lightgrey")
    ax.grid(True, linewidth=1, which="major", color="lightgrey")
    ax.axvline(x=1, linewidth=2, color='grey', linestyle="--")
    ax.legend(title="Tardigrade-proteome(s) in clusters", ncol=1, numpoints=1, loc="upper right", frameon=True, prop={'size':18})
    print "[STATUS] - Plotting %s" % (non_TA_counts_by_TA_counts_plot_f)
    f.suptitle('Counts of clusters by member-composition (hypothesis=%s)' % hypothesis, fontsize=20)
    f.savefig(non_TA_counts_by_TA_counts_plot_f, format=PLOT_FORMAT)
    plt.close()

class InputObj():
    def __init__(self, args):
        #print args
        self.cluster_f = args['--cluster_f']
        self.plot_format = args['--plot_format']
        self.domain_sources = ["GO", "IPR"] + args['--domain_sources'].split(",")

if __name__ == "__main__":
    FIGSIZE = (24,12)
    args = docopt(__doc__)
    print "[+] Start... "
    inputObj = InputObj(args)
    print "[+] Parsing domain-IDs for %s" % (",".join(inputObj.domain_sources))
    print "[+] Parsing domains in %s ..." % (inputObj.interproscan_f)

    #TA = set(['HDUJA', 'RVARI', 'ETEST', 'MTARD'])
    TA = set(['HDUJA', 'RVARI'])
    #ON = set(['EKANA', 'PCAPE', 'PSEDG'])
    ON = set([])
    AR = set(['AGAMB', 'AMELL', 'APISU', 'CLECT', 'DPOND', 'DPULE', 'ISCAP', 'NVITR', 'PHUMA', 'PXYLO', 'SINVI', 'SMARI', 'TCAST', 'TURTI'])
    #NP = set(['GORSP'])
    NP = set([])
    NE = set(['ASUUM2', 'BMALA', 'BXYLO', 'CELEG', 'MHAPL', 'PMURR', 'PPACI', 'TMURI', 'TSPIR'])
    PLOT_FORMAT = 'pdf'

    ARON = AR | ON
    NENP = NE | NP
    cluster_f = sys.argv[1]
    ALO = sys.argv[2]
    hypothesis = sys.argv[3]
    max_non_TA_count = int(sys.argv[4])
    cluster_proteome_ids_by_cluster_id = parse_clusters(cluster_f)
    non_TA_counts_by_TA_counts = calculate_results(cluster_proteome_ids_by_cluster_id, ALO)
    plot_counts(non_TA_counts_by_TA_counts)
