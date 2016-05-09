#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: prefix_fasta_based_on_field.py            [-f FASTA] [-p FIELD] [-l INT]
                                                    [-h|--help]

    Options:
        -h --help                       show this
        -f, --fasta FASTA               Fasta file
        -l, --minlen INT                Minimal length [default: 30]
        -p, --prefix FIELD              Prefix field (split on ".", start at 0)

    Description:
        - Adds prefix from filename to each sequence
        - prints sequence if longer than minimal length
"""

from __future__ import division
from docopt import docopt
import sys
from os import path

class SeqObj():
    def __init__(self, header_l, seq, prefix):
        self.header_l = header_l
        self.seq = seq
        self.prefix = prefix
        self.length = len(seq)

    def print_fixed(self):
        print ">%s.%s\n%s" % (self.prefix, self.header_l[0].replace("|",""), self.seq)

    def non_terminal_stops(self):
        non_terminal = self.seq.rstrip("*")
        return non_terminal.count("*")

def readFasta(infile):
    with open(infile) as fh:
        header, seqs = '', []
        for l in fh:
            if l[0] == '>':
                if (header):
                    yield header, ''.join(seqs)
                header, seqs = l[1:-1].split(), [] # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, ''.join(seqs)

def main(args):
    fasta_f = args['--fasta']
    prefix_idx = int(args['--prefix'])

    fasta_filename = path.basename(fasta_f)
    prefix = fasta_filename.split(".")[prefix_idx]
    for header, seq in readFasta(fasta_f):
        seqObj = SeqObj(header, seq, prefix)
        if seqObj.length >= 30:
            seqObj.print_fixed()
        #if seqObj.non_terminal_stops() > 0:
        #    seqObj.print_fixed()



if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)
