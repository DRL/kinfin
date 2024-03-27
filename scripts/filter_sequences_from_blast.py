#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage: filter_sequences_from_blast.py      -b <DIR> [-s <FILE>] [-o <DIR>] [-i <FILE>] [-e <FILE>]
                                                    [-v] [-h|--help]

    Options:
        -h --help                           show this
        -b, --blast_dir <DIR>               Directory of BLAST files to be filtered (Filenames: BlastX_Y.txt)
        -o, --output_dir <DIR>              Output directory for filtered BLAST files [default: filtered_blasts]
        -i, --include_id_f <FILE>           File containing headers of sequences to be included
        -e, --exclude_id_f <FILE>           File containing headers of sequences to be excluded
        -s, --sequence_id_f <FILE>          OrthoFinder SequenceID.txt mapping sequence headers to sequence-IDs in BLAST results.
                                                - If not provided, headers in include/exclude file are being used.
        -v, --verbose                       Verbose output
"""
from docopt import docopt
import os
import shutil
import sys
from collections import defaultdict


class DataCollection():
    def __init__(self, out_path, include_ids, exclude_ids):
        self.out_path = out_path
        self.include_ids = include_ids
        self.exclude_ids = exclude_ids
        self.status_by_blast_id = {}

    def get_status_by_blast_id(self, sequence_id_f):
        if self.exclude_ids:
            self.status_by_blast_id = defaultdict(lambda: True)
        elif self.include_ids:
            self.status_by_blast_id = defaultdict(lambda: False)
        else:
            pass
        if sequence_id_f:
            ex_in_matches_in_sequence_id_f = 0
            for line in read_file(sequence_id_f):
                try:
                    blast_id, sequence_id = line.rstrip("\n").split(": ")
                except ValueError:
                    sys.exit("[X] File %s does not seem to be in the right format." % (sequence_id_f))
                if self.exclude_ids and sequence_id in self.exclude_ids:
                    self.status_by_blast_id[blast_id] = False
                    ex_in_matches_in_sequence_id_f += 1
                if self.include_ids and sequence_id in self.include_ids:
                    self.status_by_blast_id[blast_id] = True
                    ex_in_matches_in_sequence_id_f += 1
            if ex_in_matches_in_sequence_id_f == 0:
                sys.exit("[X] No ID in exclude_f/include_f matched IDs in sequence_id_f.")
        else:
            if self.exclude_ids:
                for blast_id in self.exclude_ids:
                    self.status_by_blast_id[blast_id] = False
            if self.include_ids:
                for blast_id in self.include_ids:
                    self.status_by_blast_id[blast_id] = True

    def filter_blast_dir(self, blast_dir):
        blast_f_total = 0
        blast_f_done = 1
        for f in os.listdir(blast_dir):
            if f.startswith("Blast") and f.endswith(".txt"):
                blast_f_total += 1
        for f in os.listdir(blast_dir):
            if f.startswith("Blast") and f.endswith(".txt"):
                blast_f = os.path.join(blast_dir, f)
                print("[+] Filtering %s (%s/%s)" % (blast_f, blast_f_done, blast_f_total))
                outfile = os.path.join(self.out_path, f)
                statsfile = os.path.join(self.out_path, f + ".stats.txt")
                outstring = []
                stats_dict = {"total": 0, "excluded": 0, "included": 0}
                for line in read_file(blast_f):
                    stats_dict["total"] += 1
                    field = line.rstrip("\n").split("\t")
                    status = [self.status_by_blast_id[field[0]], self.status_by_blast_id[field[1]]]
                    if verbose_flag:
                        print(field)
                        print(status, all(status))
                    if all(status) is True:
                        outstring.append(line)
                        stats_dict["included"] += 1
                    else:
                        stats_dict["excluded"] += 1
                with open(outfile, 'w') as out_fh:
                    out_fh.write("".join(outstring))
                stats_dict["total"] = stats_dict["included"] + stats_dict["excluded"]
                if stats_dict["total"] == 0:
                    print("[-] No blast results found in %s." % (blast_f))
                    with open(statsfile, 'w') as stats_fh:
                        stats_fh.write("file=%s;total=0;FILE WAS EMPTY!.\n")
                else:
                    with open(statsfile, 'w') as stats_fh:
                        stats_fh.write("file=%s;total=%s;included=%s;included_perc=%s;excluded=%s;excluded_perc=%s\n" % (
                            outfile,
                            stats_dict["total"],
                            stats_dict["included"],
                            "{:.1%}".format(stats_dict["included"] / stats_dict["total"]),
                            stats_dict["excluded"],
                            "{:.1%}".format(stats_dict["excluded"] / stats_dict["total"])))
                blast_f_done += 1


def read_file(infile):
    if not os.path.exists(infile):
        sys.exit("[X] File provided for argument --include_id_f/--exclude_id_f does not exist.")
    else:
        with open(infile) as fh:
            for line in fh:
                yield line


def get_ids(include_id_f, exclude_id_f):
    include_ids = []
    exclude_ids = []
    if include_id_f and exclude_id_f:
        sys.exit("[X] Please provide an argument for EITHER --include_id_f or --exclude_id_f. Can't be both.")
    elif include_id_f:
        include_ids = parse_ids_f(include_id_f)
    elif exclude_id_f:
        exclude_ids = parse_ids_f(exclude_id_f)
    else:
        sys.exit("[X] Please provide an argument for --include_id_f or --exclude_id_f.")
    return set(include_ids), set(exclude_ids)


def parse_ids_f(id_f):
    ids = []
    for line in read_file(id_f):
        ids.append(line.rstrip("\n"))
    return set(ids)


def get_output_path(out_dir):
    out_path = os.path.join(os.getcwd(), out_dir)
    if os.path.exists(out_path):
        valid_decisions = set(['y', 'n'])
        decision = ""
        while decision not in valid_decisions:
            decision = input("[-] Directory %s already exists\n\tOverwrite (y/n) ? : " % (out_path))
        if decision == 'y':
            shutil.rmtree(out_path)
        else:
            sys.exit("[X] Exiting.")
    print("[+] Creating directory \n\t%s" % (out_path))
    os.mkdir(out_path)
    return out_path


if __name__ == "__main__":
    __version__ = 0.3
    args = docopt(__doc__)

    blast_dir = args['--blast_dir']
    out_dir = args['--output_dir']
    include_id_f = args['--include_id_f']
    exclude_id_f = args['--exclude_id_f']
    sequence_id_f = args['--sequence_id_f']
    verbose_flag = args['--verbose']

    # check, create out_dir and get path
    out_path = get_output_path(out_dir)

    # include/exclude
    include_ids, exclude_ids = get_ids(include_id_f, exclude_id_f)
    # Init data
    dataCollection = DataCollection(out_path, include_ids, exclude_ids)
    # parse sequence_id_f depending on include_ids/exclude_ids
    dataCollection.get_status_by_blast_id(sequence_id_f)
    # filte BLAST based on sequence_ids
    dataCollection.filter_blast_dir(blast_dir)
