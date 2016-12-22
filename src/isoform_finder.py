#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: find_isoforms.py                      -f <FASTA> -g <GFF3> -t <TYPE>
                                             [-o <STRING>]
                                              [-h|--help]

    Options:
        -h --help                           show this
        -f, --fasta <FILE>                  FASTA file
        -g, --gff3 <FILE>                   GFF3 file
        -t, --type <STRING>                 GFF3 file type (ENSEMBL, WormBase, ...)
        -o, --outprefix <STRING>            Output prefix

"""

import re
import sys
import operator
from docopt import docopt
from os.path import isfile, join, exists, realpath, dirname, basename

def read_gff(gff_f, gff_type):
    with open(gff_f) as gff_fh:
        for line in gff_fh:
            if not line.startswith("#"):
                temp = line.rstrip("\n").split("\t")
                if len(temp) == 9 and temp[1] == gff_type:
                    record_type = temp[2]
                    ninth_by_key = {}
                    for field in temp[8].split(";"):
                        if field:
                            key_value = field.split("=")
                            ninth_by_key[key_value[0]] = key_value[1]
                    yield record_type, ninth_by_key

def read_fasta_length(fasta_f):
    with open(fasta_f) as fasta_fh:
        header, seqs = '', []
        for l in fasta_fh:
            if l[0] == '>':
                if header:
                    yield header, len(''.join(seqs))
                header, seqs = l[1:-1].split()[0], [] # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        yield header, len(''.join(seqs))

def parse_fasta(fasta_f):
    proteome_id = basename(fasta_f).split(".")[0]
    proteinCollection = ProteinCollection(proteome_id)
    for protein_id, protein_length in read_fasta_length(fasta_f):
        proteinCollection.add_protein(protein_id, protein_length)
    return proteinCollection

class ProteinCollection():
    def __init__(self, proteome_id):
        self.proteome_id = proteome_id
        self.protein_ids = []
        self.protein_length_by_protein_id = {}
        self.protein_count = 0

    def add_protein(self, protein_id, protein_length):
        raw_protein_id = ".".join(protein_id.split(".")[1:])
        self.protein_ids.append(raw_protein_id)
        self.protein_length_by_protein_id[raw_protein_id] = protein_length
        self.protein_count += 1

    def get_protein_length(self, raw_protein_id):
        if raw_protein_id in self.protein_length_by_protein_id:
            protein_id = "%s.%s" % (self.proteome_id, raw_protein_id)
            return (protein_id, self.protein_length_by_protein_id[raw_protein_id])

def parse_gff(gff_f, gff_type):
    annotationCollection = AnnotationCollection()
    gene_id_by_mRNA_id = {}
    for record_type, ninth_by_key in read_gff(gff_f, gff_type):
        if gff_type == "NCBI" or gff_type == "Gnomon" or gff_type == "Ensembl_Metazoa" or gff_type == "BeetleBase" or gff_type == "NasoniaBase" or gff_type == "AphidBase" or gff_type == "ORCAE":
            #print record_type, ninth_by_key
            if record_type == 'CDS':
                mRNA_id = ninth_by_key['Parent']
                protein_id = ninth_by_key['protein_id']
                #print gene_id_by_mRNA_id
                gene_id = gene_id_by_mRNA_id[mRNA_id]
                annotationCollection.add_mRNA(gene_id, protein_id)
            elif record_type == 'mRNA':
                gene_id = ninth_by_key['Parent']
                mRNA_id = ninth_by_key['ID']
                if not gene_id in gene_id_by_mRNA_id:
                    gene_id_by_mRNA_id[mRNA_id] = gene_id
                else:
                    sys.exit("[-] Collision of mRNA_id %s." % (mRNA_id))
            else:
                pass
        else:
            if record_type == 'gene':
                gene_id = ninth_by_key['ID']
                if gene_id in annotationCollection.mRNA_ids_by_gene_id:
                    sys.exit("[-] GeneID %s duplicated in GFF3 file %s." % (gene_id, gff_f))
                else:
                    annotationCollection.add_gene(gene_id)
            elif record_type == "mRNA":
                parent_id = ninth_by_key['Parent']
                mRNA_id = None
                if gff_type == "VectorBase" or gff_type == "BeeBase" or gff_type == "I5K":
                    mRNA_id = ninth_by_key['ID'].replace("transcript:", "").replace("-R", "-P")
                elif gff_type == "JGI":
                    #mRNA_id = ninth_by_key['ID'].replace("transcript:", "").replace("T", "P") # CTELE.capitella_teleta.GCA_000328365_1.ENSEMBL.protein.gff3
                    mRNA_id = ninth_by_key['ID'].replace("transcript:", "") # DPULE.daphnia_pulex.GCA_000187875.ENSEMBLE30.protein.gff3
                elif gff_type == "TRIA":
                    mRNA_id = ninth_by_key['ID'].replace("transcript:", "")
                elif gff_type == "OIST":
                    mRNA_id = ninth_by_key['ID'].replace("transcript:", "") + ".p"
                else:
                    mRNA_id = ninth_by_key['Name']
                print "parent_id", parent_id, "mRNA_id", mRNA_id
                annotationCollection.add_mRNA(parent_id, mRNA_id)
            else:
                pass
    return annotationCollection

class AnnotationCollection():
    def __init__(self):
        self.mRNA_ids_by_gene_id = {}
        self.gene_count = 0
        self.mRNA_count = 0

    def get_longest_mRNA_id(self, gene_id):
        lengths = []
        protein_ids = []
        for mRNA_id in self.mRNA_ids_by_gene_id[gene_id]:
            if mRNA_id in proteinCollection.protein_length_by_protein_id:
                protein_ids.append(mRNA_id)
                lengths.append(proteinCollection.protein_length_by_protein_id[mRNA_id])
            else:
                truncated_mRNA_id = ".".join(mRNA_id.split(".")[0:-1])
                if truncated_mRNA_id in proteinCollection.protein_length_by_protein_id:
                    protein_ids.append(truncated_mRNA_id)
                    lengths.append(proteinCollection.protein_length_by_protein_id[truncated_mRNA_id])
                #else:
                #    protein_ids.append(mRNA_id)
                #    lengths.append(0)
        if lengths:
            mRNA_index, protein_length = max(enumerate(lengths), key=operator.itemgetter(1))
            return protein_ids[mRNA_index], protein_length, list(set(protein_ids))
        else:
            return None, None, None

    def add_gene(self, gene_id):
        self.mRNA_ids_by_gene_id[gene_id] = set()
        self.gene_count += 1

    def add_mRNA(self, gene_id, mRNA_id):
        if not gene_id in self.mRNA_ids_by_gene_id:
            self.mRNA_ids_by_gene_id[gene_id] = set()
            self.gene_count += 1
        self.mRNA_ids_by_gene_id[gene_id].add(mRNA_id)
        self.mRNA_count += 1

def generate_output(out_prefix):
    out_longest_isoforms_f = ""
    out_non_longest_isoforms_f = ""
    out_mRNAs_not_in_fasta_f = ""
    if out_prefix:
        out_non_longest_isoforms_f = "%s.%s.non_longest_isoforms.txt" % (out_prefix, proteinCollection.proteome_id)
        out_longest_isoforms_f = "%s.%s.longest_isoforms.txt" % (out_prefix, proteinCollection.proteome_id)
        out_mRNAs_not_in_fasta_f = "%s.%s.mRNAs_not_in_fasta.txt" % (out_prefix, proteinCollection.proteome_id)
        out_proteins_not_in_gff_f = "%s.%s.proteins_not_in_gff.txt" % (out_prefix, proteinCollection.proteome_id)
    else:
        out_non_longest_isoforms_f = "%s.non_longest_isoforms.txt" % (proteinCollection.proteome_id)
        out_longest_isoforms_f = "%s.longest_isoforms.txt" % (proteinCollection.proteome_id)
        out_mRNAs_not_in_fasta_f = "%s.mRNAs_not_in_fasta.txt" % (proteinCollection.proteome_id)
        out_proteins_not_in_gff_f = "%s.proteins_not_in_gff.txt" % (proteinCollection.proteome_id)

    longest_isoforms = []
    non_longest_isoforms = []
    mRNAs_not_in_fasta = []
    fasta_proteins_seen = set()

    for gene_id, mRNA_ids in annotationCollection.mRNA_ids_by_gene_id.items():
        if len(mRNA_ids):
            longest_protein_id, longest_protein_length, protein_ids = annotationCollection.get_longest_mRNA_id(gene_id)
            #print longest_protein_id, longest_protein_length, protein_ids
            if longest_protein_length:
                longest_isoforms.append("%s.%s" % (proteinCollection.proteome_id, longest_protein_id))
                for protein_id in protein_ids:
                    fasta_proteins_seen.add(protein_id)
                    if not protein_id == longest_protein_id:
                        non_longest_isoforms.append("%s.%s" % (proteinCollection.proteome_id, protein_id))
            else:
                for mRNA_id in mRNA_ids:
                    mRNAs_not_in_fasta.append("%s.%s" % (proteinCollection.proteome_id, mRNA_id))


    proteins_not_in_gff = ["%s.%s" % (proteinCollection.proteome_id, protein_id) for protein_id in list(set(proteinCollection.protein_length_by_protein_id.keys()).difference(fasta_proteins_seen))]

    print "[+] %s proteinIDs in FASTA file" % (proteinCollection.protein_count)
    if longest_isoforms:
        print "[+] Writing %s proteinIDs to %s ..." % (len(longest_isoforms), out_longest_isoforms_f)
        with open(out_longest_isoforms_f, 'w') as out_longest_isoforms_fh:
            out_longest_isoforms_fh.write("\n".join(longest_isoforms) + "\n")
    if non_longest_isoforms:
        print "[+] Writing %s proteinIDs to %s ..." % (len(non_longest_isoforms), out_non_longest_isoforms_f)
        with open(out_non_longest_isoforms_f, 'w') as out_non_longest_isoforms_fh:
            out_non_longest_isoforms_fh.write("\n".join(non_longest_isoforms) + "\n")
    if mRNAs_not_in_fasta:
        print "[+] Writing %s proteinIDs to %s ..." % (len(mRNAs_not_in_fasta), out_mRNAs_not_in_fasta_f)
        with open(out_mRNAs_not_in_fasta_f, 'w') as out_mRNAs_not_in_fastah_f:
            out_mRNAs_not_in_fastah_f.write("\n".join(mRNAs_not_in_fasta) + "\n")
    if proteins_not_in_gff:
        print "[+] Writing %s proteinIDs to %s ..." % (len(proteins_not_in_gff), out_proteins_not_in_gff_f)
        with open(out_proteins_not_in_gff_f, 'w') as out_proteins_not_in_gff_fh:
            out_proteins_not_in_gff_fh.write("\n".join(proteins_not_in_gff) + "\n")

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)

    fasta_f = args['--fasta']
    gff_f = args['--gff3']
    gff_type = args['--type']
    out_prefix = args['--outprefix']

    proteinCollection = None
    print "[+] Start ..."
    if fasta_f and gff_f:
        print "[+] Parsing FASTA %s ..." % fasta_f
        proteinCollection = parse_fasta(fasta_f)
        print "[+] Parsing GFF3 %s ..." % gff_f
        annotationCollection = parse_gff(gff_f, gff_type)
        generate_output(out_prefix)
    else:
        print "[-] No input files given."
        sys.exit(__doc__.strip())

