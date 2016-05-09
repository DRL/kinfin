#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import os
import gzip

def parse_sequence_ids_to_exclude(sequence_ids_f):
    exclude_set = set()
    with open(sequence_ids_f) as fh:
        for l in fh:
            sequence_id = l.rstrip("\n").split(":")[0]
            exclude_set.add(sequence_id)
    return exclude_set

def filter_blast_dir(blast_f):

    if blast_f.endswith(".gz"):
        infile = blast_f
        temp = os.path.basename(blast_f)
        print "IN: ", infile
        outfile = os.path.join(output_dir, temp)
        print "OUT: ", outfile
        output = ''
        with gzip.open(infile) as in_fh:
            for line in in_fh:
                temp = line.rstrip("\n").split()
                if temp[0] in exclude_set or temp[1] in exclude_set:
                    pass
                else:
                    output += line
            if (output):
                with gzip.open(outfile, "w") as out_fh:
                    out_fh.write(output)

if __name__ == "__main__":
    try:
        blast_file = sys.argv[1]
        sequence_ids_f = sys.argv[2]
        output_dir = sys.argv[3]
    except:
        sys.exit("script.py BLASTFILE SEQUENCEIDSTOEXCLUDE OUTPUTDIR")

    exclude_set = parse_sequence_ids_to_exclude(sequence_ids_f)
    filter_blast_dir(blast_file)

