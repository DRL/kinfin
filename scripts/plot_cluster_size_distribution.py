#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: plot_cluster_size_distribution.py    -i <FILE> [-o <STRING>] [-p <STRING>] [-c <STRING>] [-x <INT>]
                                            [-h|--help]

    Options:
        -h --help                           show this
        -i, --infile <FILE>                 cluster_counts_by_taxon.txt from KinFin analysis
        -o, --out_prefix <STRING>           Outprefix [default: cluster_size_distribution]
        -c, --colormap <STRING>             Matplotlib colormap name [default: Paired]
        -x, --xlim <INT>                    xlim for "loglin", "logbar" and "barperc" plot [default: 200]
        -p, --plot_fmt <STRING>             Plot format [default: png]

"""
from __future__ import division
import re
import sys
import powerlaw
from docopt import docopt
from collections import OrderedDict
from os.path import isfile, join, exists, realpath, dirname, basename
import matplotlib as mat
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
mat.use("agg")
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def progress(iteration, steps, max_value):
    if int(iteration) == int(max_value):
        sys.stdout.write('\r')
        print "\t[PROGRESS] \t- %d%%" % (100)
    elif int(iteration) % int(steps+1) == 0:
        sys.stdout.write('\r')
        print "\t [PROGRESS] \t- %d%%" % (float(int(iteration)/int(max_value))*100),
        sys.stdout.flush()
    else:
        pass

def read_file(infile):
    if not infile or not exists(infile):
        sys.exit("[ERROR] - File '%s' does not exist." % (infile))
    if infile.endswith(".gz"):
        import gzip
        with gzip.open(infile) as fh:
            for line in fh:
                yield line.rstrip("\n")
    else:
        with open(infile) as fh:
            for line in fh:
                yield line.rstrip("\n")

class DataObj():
    def __init__(self, out_prefix, plot_fmt, cmap, xlim):
        self.out_prefix = out_prefix
        self.plot_fmt = plot_fmt
        self.cmap = cmap
        self.xlim = xlim
        self.counts_by_cluster_size_by_proteome_counts = {}
        self.counts_by_cluster_size = {}
        self.cluster_sizes = []
        self.cluster_size_max = 0
        self.cluster_count_max = 0
        self.proteome_count_max = 0
        self.count_of_cluster_sizes = 0
        self.proteome_counts_by_colour = OrderedDict()

    def parse_data(self, infile):
        for line in read_file(infile):
            if not line.startswith("#") and not line.startswith("ID"):
                temp = line.split("\t")
                cluster_size = sum([int(count) for count in temp[1:]])
                proteome_count = sum([1 for count in temp[1:] if not count == '0'])
                self.add_cluster(proteome_count, cluster_size)
        self.homogenise_counts()

    def add_cluster(self, proteome_count, cluster_size):
        if not proteome_count in self.counts_by_cluster_size_by_proteome_counts:
            self.counts_by_cluster_size_by_proteome_counts[proteome_count] = {}
        self.counts_by_cluster_size_by_proteome_counts[proteome_count][cluster_size] = self.counts_by_cluster_size_by_proteome_counts[proteome_count].get(cluster_size, 0) + 1
        self.counts_by_cluster_size[cluster_size] = self.counts_by_cluster_size.get(cluster_size, 0) + 1

    def homogenise_counts(self):
        for proteome_count in self.counts_by_cluster_size_by_proteome_counts:
            for cluster_size in self.counts_by_cluster_size:
                if not cluster_size in self.counts_by_cluster_size_by_proteome_counts[proteome_count]:
                    self.counts_by_cluster_size_by_proteome_counts[proteome_count][cluster_size] = 0
        self.proteome_count_max = max(list(self.counts_by_cluster_size_by_proteome_counts.keys()))
        self.count_of_cluster_sizes = len(self.counts_by_cluster_size)
        self.cluster_sizes = sorted(list(self.counts_by_cluster_size.keys()))
        self.cluster_size_max = max(self.cluster_sizes)
        self.cluster_count_max = max(self.counts_by_cluster_size.values())

    def yield_counts_by_proteome_count(self, count_type):
        self.proteome_counts_by_colour = OrderedDict()
        cmap = mat.cm.get_cmap(self.cmap)
        y_bottom = [0] * self.count_of_cluster_sizes
        x = self.cluster_sizes
        for proteome_count in sorted(self.counts_by_cluster_size_by_proteome_counts):
            colour = cmap(1.*proteome_count/self.proteome_count_max)
            if not colour in self.proteome_counts_by_colour:
                self.proteome_counts_by_colour[colour] = []
            self.proteome_counts_by_colour[colour].append(proteome_count)
            if count_type == "absolute":
                y = [self.counts_by_cluster_size_by_proteome_counts[proteome_count][cluster_size] + y_b for cluster_size, y_b in zip(sorted(self.counts_by_cluster_size_by_proteome_counts[proteome_count]), y_bottom)]
                #filled = self.interpolate_gaps(y, 1)
                yield proteome_count, x, y, y_bottom, colour
                #y_bottom = [_y if _y >= 1 else 0.01 for _y in y]
                y_bottom = y
            elif count_type == "relative":
                y = [self.counts_by_cluster_size_by_proteome_counts[proteome_count][cluster_size] for cluster_size in sorted(self.counts_by_cluster_size_by_proteome_counts[proteome_count])]
                yield proteome_count, x, y, y_bottom, colour
                y_bottom = [y_ + y_b for y_, y_b in zip(y, y_bottom)]
            else:
                pass

    def get_cluster_size_and_count(self):
        x, y = [], []
        for cluster_size, count in sorted(self.counts_by_cluster_size.items()):
            x.append(cluster_size)
            y.append(count)
        return x, y

    def plot_cluster_sizes(self, plot_type):
        f, ax = plt.subplots(figsize=(12,6))
        out_f = ''
        ax.set_facecolor('white')
        print "[+] Plotting \"%s\" ..." % (plot_type)
        if plot_type in PLOTS:
            if plot_type == "loglog" or plot_type == "loglin" or plot_type == "loglogpowerlaw":
                for proteome_count, x, y, y_bottom, colour in self.yield_counts_by_proteome_count('absolute'):
                    ax.plot(x, y, c='None')
                    mask = np.where(y >= 1, y, "nan")
                    ax.fill_between(x, y, y_bottom, facecolor=colour, interpolate = True, where = mask)
                    progress(proteome_count, 1, self.proteome_count_max)
                ax.plot
            elif plot_type == "logbar":
                ax.set_yscale('log')
                for proteome_count, x, y, y_bottom, colour in self.yield_counts_by_proteome_count('relative'):
                    ax.bar(x, y, bottom=y_bottom, color=colour)
                    progress(proteome_count, 1, self.proteome_count_max)
            elif plot_type == "barperc":
                for proteome_count, x, y, y_bottom, colour in self.yield_counts_by_proteome_count('absolute'):
                    y_perc = [(y_ / self.counts_by_cluster_size[cluster_size]) * 100 for y_, cluster_size in zip(y, sorted(self.counts_by_cluster_size))]
                    y_bottom_perc = [(y_b / self.counts_by_cluster_size[cluster_size]) * 100 for y_b, cluster_size in zip(y_bottom, sorted(self.counts_by_cluster_size))]
                    ax.bar(x, y_perc, bottom=y_bottom_perc, color=colour)
                    progress(proteome_count, 1, self.proteome_count_max)
            else:
                pass
            if plot_type == "loglog":
                out_f = "%s.loglog.%s" % (self.out_prefix, self.plot_fmt)
                ax.set_xscale('log')
                plt.gca().set_ylim(bottom=0.8, top=self.cluster_count_max * 2)
                plt.gca().set_xlim(left=0.8, right=self.cluster_size_max * 2)
            elif plot_type == "loglin":
                out_f = "%s.loglin.%s" % (self.out_prefix, self.plot_fmt)
                plt.gca().set_ylim(bottom=0.8, top=self.cluster_count_max * 2)
                plt.gca().set_xlim(left=-2, right=self.xlim)
            elif plot_type == "logbar":
                out_f = "%s.logbar.%s" % (self.out_prefix, self.plot_fmt)
                plt.gca().set_ylim(bottom=0.8, top=self.cluster_count_max * 2)
                plt.gca().set_xlim(left=-2, right=self.xlim)
            elif plot_type == "barperc":
                out_f = "%s.barperc.%s" % (self.out_prefix, self.plot_fmt)
                plt.gca().set_ylim(bottom=0.8, top=100)
                plt.gca().set_xlim(left=-2, right=self.xlim)
            else:
                if plot_type == "powerlaw" or plot_type == "loglogpowerlaw":
                    x, y = self.get_cluster_size_and_count()
                    x_log_array = np.log10(np.array(x))
                    y_log_array = np.log10(np.array(y))
                    x_array = np.array(x)
                    y_array = np.array(y)
                    fit = np.polyfit(x_log_array, y_log_array, 1)
                    powerlaw = lambda x, amp, index: amp * (x**index)
                    powerlaw_y = powerlaw(x_array, fit[1]*y[1], fit[0])
                    ax.plot(x_array, powerlaw_y, '--', c='grey', alpha=0.8, ms=5, label="Power-Law")
                    #################
                    ax.plot(x_array, y_array, 'o', c='black', ms=2, alpha=0.5, label="Clustering")
                    if plot_type == "powerlaw":
                        out_f = "%s.powerlaw.%s" % (self.out_prefix, self.plot_fmt)
                    else:
                        out_f = "%s.loglogpowerlaw.%s" % (self.out_prefix, self.plot_fmt)
                    #ax.legend()
                    plt.gca().set_ylim(bottom=0.9, top=self.cluster_count_max * 2)
                    plt.gca().set_xlim(left=0.8, right=self.cluster_size_max * 2)
                    ax.set_yscale('symlog', linthreshy=1.0)
                    ax.set_xscale('log')
            legend_handles = []
            for colour, proteomes in self.proteome_counts_by_colour.items():
                if len(proteomes) > 1:
                    legend_handles.append(mpatches.Patch(label="%s-%s" % (proteomes[0], proteomes[-1]), color=colour))
                else:
                    legend_handles.append(mpatches.Patch(label="%s" % (proteomes[0]), color=colour))
            legend = ax.legend(handles=legend_handles, ncol=2, loc='best', numpoints=1, frameon=True, title="Number of proteomes in cluster")
            legend.get_frame().set_facecolor('white')
            plt.margins(1)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax.grid(True, linewidth=0.5, which="minor", color="lightgrey")
            ax.grid(True, linewidth=1, which="major", color="lightgrey")
            if plot_type == "barperc":
                ax.set_xlabel('Cluster size')
                ax.set_ylabel('Number of proteomes in cluster (%)')
            else:
                ax.set_xlabel('Cluster size')
                ax.set_ylabel('Count')
            f.tight_layout()
            f.savefig(out_f, format=self.plot_fmt)
            plt.close()


if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    input_f = args['--infile']
    out_prefix = args['--out_prefix']
    cmap = args['--colormap']
    plot_fmt = args['--plot_fmt']
    xlim = int(args['--xlim'])

    PLOTS = set(['loglog', 'loglin', 'logbar', 'barperc', 'loglogpowerlaw'])
    print "[+] Start ..."
    dataObj = DataObj(out_prefix, plot_fmt, cmap, xlim)
    dataObj.parse_data(input_f)
    #dataObj.plot_cluster_sizes('loglog')
    #dataObj.plot_cluster_sizes('loglin')
    #dataObj.plot_cluster_sizes('logbar')
    #dataObj.plot_cluster_sizes('barperc')
    #dataObj.plot_cluster_sizes('powerlaw')
    dataObj.plot_cluster_sizes('loglogpowerlaw')
