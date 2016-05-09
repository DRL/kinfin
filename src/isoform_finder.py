#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import os

'''
A : longest, shortest isoform
- where do they live in the clusters?
- what is the influence of analysing a subset of these?

B: stats file with isoforms per gene
- gene, length, isoform_count, isoform_list
'''

class DataObj():
    def __init__(self):
        self.order_of_fasta_files = []
        self.fasta_files_by_id = {}

    def add_annotation(self, annotationObj):
        fasta_f = annotationObj.fasta_f
        self.order_of_fasta_files.append(fasta_f)
        self.fasta_files_by_id[fasta_f] = annotationObj

    def output(self):
        data_stats_f = "isoform_stats.txt"
        with open(data_stats_f, "w") as fh:
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (\
                    "species_id", \
                    "species_name", \
                    "fasta_f", \
                    "count_genes", \
                    "count_sequences", \
                    "count_isoforms" \
                    ))
            for fasta_f in self.order_of_fasta_files:
                annotationObj = self.fasta_files_by_id[fasta_f]

                annotationObj.output("longest")

                species_id = annotationObj.species_id
                species_name = annotationObj.species_name
                count_genes = annotationObj.count_genes
                count_sequences = annotationObj.count_sequences
                count_isoforms = annotationObj.count_isoforms
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (\
                    species_id, \
                    species_name, \
                    fasta_f, \
                    count_genes, \
                    count_sequences, \
                    count_isoforms \
                    ))

class AnnotationObj():
    def __init__(self, fasta_f):
        self.order_of_genes = []
        self.genes_by_name = {}
        self.count_genes = 0
        self.count_sequences = 0
        self.count_isoforms = 0
        self.fasta_f = fasta_f
        self.species_id = fasta_f.split(".")[0]
        self.species_name = fasta_f.split(".")[1]

    def output(self, criterion):

        fasta_stats_f = self.fasta_f + ".stats.txt"
        fasta_included_f = self.fasta_f + "." + criterion + ".faa"
        fasta_excluded_f = self.fasta_f + ".not_" + criterion + ".faa"

        fasta_stats_fh = open(fasta_stats_f, "w")

        fasta_excluded_string = ''
        fasta_included_string = ''

        for geneObj in self.yield_geneObjs():

            seqObj = ''
            if geneObj.isoforms() >= 1:
                if criterion == "longest":
                    seqObj = geneObj.get_longest_isoform()
                elif criterion == "shortest":
                    seqObj = geneObj.get_shortest_isoform()
                else:
                    sys.exit("[ERROR2]")
            else:
                seqObj = geneObj.get_sequence()
            seqObj.criterion = 1
            fasta_included_string += ">%s\n%s\n" % (seqObj.header, seqObj.sequence)

            geneObj_stats = geneObj.get_stats()
            fasta_stats_fh.write(geneObj_stats)

            for seqObj in geneObj.sequences:
                if seqObj.criterion == 0:
                    fasta_excluded_string += ">%s\n%s\n" % (seqObj.header, seqObj.sequence)

        if (fasta_included_string):
            fasta_included_fh = open(fasta_included_f, "w")
            fasta_included_fh.write(fasta_included_string)
            fasta_included_fh.close()

        if (fasta_excluded_string):
            fasta_excluded_fh = open(fasta_excluded_f, "w")
            fasta_excluded_fh.write(fasta_excluded_string)
            fasta_excluded_fh.close()

        fasta_stats_fh.close()


    def add_sequence(self, seqObj):
        gene_name = seqObj.gene
        if not gene_name in self.genes_by_name:

            self.count_genes += 1
            self.order_of_genes.append(gene_name)
            geneObj = GeneObj(gene_name)

            self.genes_by_name[gene_name] = GeneObj(gene_name)
        else:
            self.count_isoforms += 1
        self.count_sequences += 1
        self.genes_by_name[gene_name].add_seq(seqObj)

    def yield_geneObjs(self):
        for gene_name in self.order_of_genes:
            yield self.genes_by_name[gene_name]


class GeneObj():
    def __init__(self, name):
        self.name = name
        self.sequences = []

    def get_stats(self):
        string = ''
        for seqObj in self.sequences:
            this_header = seqObj.short_header
            isoforms = ",".join([x.short_header for x in self.sequences if not x.short_header == this_header])
            string += "%s\t%s\t%s\t%s\t%s\t" % (\
                    seqObj.short_header, \
                    seqObj.gene, \
                    seqObj.length, \
                    seqObj.criterion, \
                    self.isoforms())
            if (isoforms):
                string += "%s\n" % isoforms
            else:
                string += "\n"
        return string


    def add_seq(self, seqObj):
        self.sequences.append(seqObj)

    def get_longest_isoform(self):
        longest = max([x.length for x in self.sequences])
        for seqObj in self.sequences:
            if seqObj.length == longest:
                return seqObj

    def get_shortest_isoform(self):
        shortest = min([x.length for x in self.sequences])
        for seqObj in self.sequences:
            if seqObj.length == shortest:
                return seqObj

    def get_sequence(self):
        return self.sequences[0]

    def isoforms(self):
        return len(self.sequences) - 1

class SeqObj():
    def __init__(self, gene, header, sequence):
        self.header = header
        self.short_header = header.split()[0]
        self.sequence = sequence
        self.length = len(sequence)
        self.gene = gene
        self.criterion = 0

def parse_fasta(fasta_f):
    header, sequence = '', ''
    regex = REGEX_D[fasta_f]

    annotationObj = AnnotationObj(fasta_f)

    with open(fasta_f) as fh:
        i = 0
        for l in fh:
            line = l.rstrip("\n")
            if line.startswith(">"):
                i += 1
                # deal with previous header
                if (header):
                    gene = get_gene_name(header, regex)
                    seqObj = SeqObj(gene, header, sequence)
                    annotationObj.add_sequence(seqObj)
                    #if i == 100:
                    #    break
                header = line.lstrip(">")
                sequence = ''
            else:
                sequence += line
        gene = get_gene_name(header, regex)
        seqObj = SeqObj(gene, header, sequence)
        annotationObj.add_sequence(seqObj)
    return annotationObj

def get_gene_name(header, regex):
    gene = ''
    regex = regex.strip("\"")
    #print "regex: ", regex
    #print "header: ", header
    if regex == "NONE":
        gene = header
    elif regex == "column1/augustus":
        gene = ".".join(header.split(".")[0:-1])
    elif regex == "column1/split('|')[0]/split('_')[2]/complete":
        gene = "_".join(header.split("|")[0].split("_")[0:-1]) + "_seq"
    elif regex == "column1/split('|')[0:1]/split('_')[2]/complete":
        temp = "|".join(header.split("|")[0:2])
        gene = "_".join(temp.split("_")[0:-1]) + "_seq"
    elif regex == "column3":
        gene = header.split()[2]
    elif regex == "gene:":
        match = re.search(r"gene:(\w+)", header)
        gene = match.group(1)
    elif regex == "gene_id=":
        match = re.search(r"gene_id=(\S+)", header)
        gene = match.group(1)
    elif regex == "parent=/split(',')":
        match = re.search(r"parent=(\w+),", header)
        gene = match.group(1)
    else:
        sys.exit("[ERROR6]")
    #print "gene: ", gene
    if not (gene) or (len(gene) < 5):
        if not gene.startswith("g"):
            sys.exit(header)

    return gene

def parse_regex(regex_f):
    regex_d = {}
    order_of_species = []
    with open(regex_f) as fh:
        for l in fh:
            line = l.rstrip("\n")
            if (line) and not line.startswith("#"):
                regex, fasta_f = line.split()
                regex_d[fasta_f] = regex
                order_of_species.append(fasta_f)
    return regex_d, order_of_species

def parse_fasta_dir(fasta_dir):
    dataObj = DataObj()
    fasta_f_in_dir = os.listdir(fasta_dir)
    for fasta_f in order_of_species:
        if fasta_f in fasta_f_in_dir:
            print fasta_f
            annotationObj = parse_fasta(fasta_f)
            dataObj.add_annotation(annotationObj)
    return dataObj

if __name__ == "__main__":
    try:
        fasta_dir = sys.argv[1]
        regex_f = sys.argv[2]
    except:
        sys.exit("script.py FASTADIR REGEXFILE")

    REGEX_D, order_of_species = parse_regex(regex_f)
    dataObj = parse_fasta_dir(fasta_dir)
    dataObj.output()
