#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import os
import gzip

def parse_sequence_ids_to_exclude(sequence_ids_f):
    exclude_set = []
    with open(sequence_ids_f) as fh:
        for l in fh:
            sequence_id = l.rstrip("\n").split(":")[0]
            exclude_set.append(sequence_id)
    return set(exclude_set)

def read_fasta(infile):
    with open(infile) as fh:
        header, seqs = '', []
        for l in fh:
            if l[0] == '>':
                if (header):
                    yield header, ''.join(seqs)
                header, seqs = l[1:-1], [] # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, ''.join(seqs)

def filter_fasta_f(fasta_f, exclude_set):
    if fasta_f.endswith(".fa"):
        infile = fasta_f
        temp = os.path.basename(fasta_f)
        print "IN: ", infile
        outfile = os.path.join(output_dir, temp)
        print "OUT: ", outfile
        output = ''
        for header, sequence in read_fasta(infile):
            if not header in exclude_set:
                output += ">%s\n%s\n" % (header, sequence)
        if (output):
            with open(outfile, "w") as out_fh:
                out_fh.write(output)

if __name__ == "__main__":
    try:
        fasta_file = sys.argv[1]
        sequence_ids_f = sys.argv[2]
        output_dir = sys.argv[3]
    except:
        sys.exit("script.py BLASTFILE SEQUENCEIDSTOEXCLUDE OUTPUTDIR")

    exclude_set = parse_sequence_ids_to_exclude(sequence_ids_f)
    filter_fasta_f(fasta_file, exclude_set)

