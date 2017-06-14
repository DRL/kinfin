#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: filter_fastas_before_clustering.py            -f FILE [-p INT] [-l INT]
                                                        [-h|--help]

    Options:
        -h --help                       show this
        -f, --fasta FILE                Fasta dir containing FASTA files
                                            - must have taxon prefix (separated by ".", beginning of filename)
                                            - must end in *.faa
        -l, --minlen INT                Minimal length [default: 30]
        -p, --prefix INT                Prefix field (split on ".", start at 0) [default: 0]

    Description:
        - Adds prefix from filename to each sequence
        - prints sequence if longer than minimal length and contains no non-terminal stops (*)
"""

from __future__ import division
from docopt import docopt
import sys
from os import path

class SeqObj():
    def __init__(self, header_l, seq, prefix):
        print header_l
        self.header = header_l[0].replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_") # orthofinder replaces chars
        print self.header
        self.seq = seq
        self.prefix = prefix
        self.length = len(seq)

    def print_fixed(self):
        print ">%s.%s\n%s" % (self.prefix, self.header, self.seq)

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
                header, seqs = l[1:-1].split(" "), [] # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, ''.join(seqs)

def main(args):
    fasta_f = args['--fasta']
    prefix_idx = int(args['--prefix'])
    min_len = int(args['--minlen'])

    fasta_filename = path.basename(fasta_f)
    prefix = fasta_filename.split(".")[prefix_idx]
    stats = {}
    stats['raw'] = 0
    stats['stop_fail'] = 0
    stats['length_fail'] = 0
    stats['pass'] = 0
    for header, seq in readFasta(fasta_f):
        stats['raw'] = stats.get('raw', 0) + 1
        seqObj = SeqObj(header, seq, prefix)
        stop_pass = 0
        length_pass = 0
        if seqObj.non_terminal_stops() > 0:
            stats['stop_fail'] += 1
        else:
            stop_pass = 1
        if seqObj.length < min_len:
            stats['length_fail'] += 1
        else:
            length_pass = 1
        if length_pass and stop_pass:
            seqObj.print_fixed()
            stats['pass'] += 1

    with open(fasta_f + ".filter_stats.txt", "w") as out_fh:
        out_fh.write("filename=%s; raw=%s (%s); stop_fail=%s (%s); length_fail=%s (%s); pass=%s (%s)\n" % (
            fasta_f,
            stats['raw'], "{:.2%}".format(stats['raw']/stats['raw']),
            stats['stop_fail'], "{:.2%}".format(stats['stop_fail']/stats['raw']),
            stats['length_fail'], "{:.2%}".format(stats['length_fail']/stats['raw']),
            stats['pass'], "{:.2%}".format(stats['pass']/stats['raw'])
            )
        )

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)
