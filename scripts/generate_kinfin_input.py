#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Usage:
        generate_kinfin_input.py -f <DIR> [-o <STRING>] [orthofinder] [config]


Options:
        -h --help                       show this
        -f, --fasta_dir <DIR>           Directory containing FASTA files ("*.fasta", "*.fas", "*.faa")
        -o, --outprefix <STRING>        Output prefix
        orthofinder                     Formats sequence header as would OrthoFinder
                                           ([:,()] are replaced with "_")
        config                          Only generates KinFin config file


"""

from __future__ import division
from docopt import docopt
import sys
import os

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    print "[+] Parsing file %s ..." % (infile)
    with open(infile) as fh:
        for line in fh:
            yield line.rstrip("\n")


def write_file(out_f, outprefix, header, lines):
    if outprefix:
        if outprefix.endswith("/"):
            if not os.path.exists(outprefix):
                os.mkdir(outprefix)
            out_f = "%s" % os.path.join(outprefix, out_f)
        else:
            out_f = "%s.%s" % (outprefix, out_f)
    print "[+] \t Writing file %s ..." % (out_f)
    with open(out_f, 'w') as out_fh:
        if header:
            out_fh.write("%s\n" % (header))
        out_fh.write("%s\n" % "\n".join(lines))


class DataCollection():
    def __init__(self, fasta_dir, outprefix, orthofinder_flag, config_flag):
        self.fasta_dir = fasta_dir
        self.outprefix = outprefix
        self.orthofinder_flag = orthofinder_flag
        self.config_flag = config_flag
        self.species_id_lines = []
        self.sequence_id_lines = []
        self.config_lines = []
        self.parse_fasta_dir()
        self.write_files()

    def parse_fasta_dir(self):
        species_idx = 0
        allowed_extensions = set(["fasta", "fas", "faa"])
        fasta_dir = self.fasta_dir
        for f in os.listdir(fasta_dir):
            if f.split(".")[-1] in allowed_extensions:
                field = f.split(".")
                self.config_lines.append("%s,%s" % (species_idx, ".".join(field[0:-1])))
                if not self.config_flag:
                    self.species_id_lines.append("%s: %s" % (species_idx, f))
                    fasta_f = os.path.join(fasta_dir, f)
                    self.parse_fasta_f(fasta_f, species_idx)
                species_idx += 1
        if not self.config_flag:
            print "[+] [Summary] %s fasta files (containing %s sequences) parsed." % (len(self.species_id_lines), len(self.sequence_id_lines))

    def parse_fasta_f(self, fasta_f, species_idx):
        seq_count = 0
        for line in read_file(fasta_f):
            if line and line[0] == '>':
                header = line[1:].split()[0]
                if self.orthofinder_flag:
                    header = header.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_")
                self.sequence_id_lines.append("%s_%s: %s" % (species_idx, seq_count, header))
                seq_count += 1
        print "[+] \t %s sequences parsed" % (seq_count)

    def write_files(self):
        print "[+] Writing output ..."
        if not self.config_flag:
            write_file('SpeciesIDs.txt', self.outprefix, None, self.species_id_lines)
            write_file('SequenceIDs.txt', self.outprefix, None, self.sequence_id_lines)
        write_file('config.txt', self.outprefix, "#IDX,TAXON", self.config_lines)


if __name__ == "__main__":
    __version__ = 0.1

    args = docopt(__doc__)
    fasta_dir = args['--fasta_dir']
    outprefix = args['--outprefix']
    orthofinder_flag = args['orthofinder']
    config_flag = args['config']
    dataCollection = DataCollection(fasta_dir, outprefix, orthofinder_flag, config_flag)
