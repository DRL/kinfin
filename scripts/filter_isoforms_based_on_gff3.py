#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: find_isoforms.py                      -f <FILE> [-g <FILE>] -t <TYPE>
                                                [--fs <INT>] [--fe <INT>] [--fdel <STRING>]
                                                [-m <FILE>] [-o <STRING>]
                                                [-h|--help]

    Options:
        -h --help                                   show this
        -f, --fasta <FILE>                          FASTA file
        --fs <INT>                                  first field of fasta header to use for GFF [default: 1]
        --fe <INT>                                  last field of fasta header to use for GFF
        --fdel <STRING>                             delimiter for fields in fasta header to use for GFF (no quotes) [default: .]
        -g, --gff3 <FILE>                           GFF3 file
        -m, --mapping <FILE>                        Mapping file of IDs (Format: FASTA-ID,GFF-ID)
        -t, --type <STRING>                         GFF3 file type (ENSEMBL, WormBase, ...) [default: None]
        -o, --outprefix <STRING>                    Output prefix

"""
from __future__ import division
import re
import sys
import operator
from docopt import docopt
from os.path import isfile, join, exists, realpath, dirname, basename

def read_gff(gff_f):
    with open(gff_f) as gff_fh:
        for line in gff_fh:
            if not line.startswith("#"):
                temp = line.rstrip("\n").split("\t")
                if len(temp) == 9:
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

def parse_mapping_dict(mapping_f):
    mapping_dict = {}
    with open(mapping_f) as mapping_fh:
        for l in mapping_fh:
            fasta_id, gff_id = l[0:-1].split(",")
            mapping_dict[fasta_id] = gff_id
    return mapping_dict

def parse_custom_fasta(fasta_f, fs, fe, fdel):
    proteome_id = basename(fasta_f).split(".")[0]
    proteinCollection = ProteinCollection(proteome_id)
    data = {}
    i = 0
    print "## Debug ..."
    for protein_id, protein_length in read_fasta_length(fasta_f):
        proteinCollection.protein_count += 1
        gene_id = ''
        if not fe:
            gene_id = eval(fdel).join(protein_id.split(eval(fdel))[fs:]) # this is used for gff
        else:
            gene_id = eval(fdel).join(protein_id.split(eval(fdel))[fs:fe]) # this is used for gff
        if not gene_id in data:
            data[gene_id] = {}
        if i <= DEBUG_ENTRIES:
            print "# FASTA-ID=%s, gene-ID=%s, LEN=%s" % (protein_id, gene_id, protein_length)
        elif i == DEBUG_ENTRIES + 1:
            print "# ..."
        else:
            pass
        i += 1
        data[gene_id][protein_id] = protein_length
    out_non_longest_isoforms_f = "%s.non_longest_isoforms.txt" % (proteinCollection.proteome_id)
    out_longest_isoforms_f = "%s.longest_isoforms.txt" % (proteinCollection.proteome_id)
    out_stats_f = "%s.stats.txt" % (proteinCollection.proteome_id)
    out_all_isoforms_by_gene_f = "%s.isoforms_by_gene_id.txt" % (proteinCollection.proteome_id)
    if out_prefix:
        out_non_longest_isoforms_f = "%s.%s.non_longest_isoforms.txt" % (out_prefix, proteinCollection.proteome_id)
        out_longest_isoforms_f = "%s.%s.longest_isoforms.txt" % (out_prefix, proteinCollection.proteome_id)
        out_stats_f = "%s.%s.stats.txt" % (out_prefix, proteinCollection.proteome_id)
        out_all_isoforms_by_gene_f = "%s.%s.isoforms_by_gene_id.txt" % (out_prefix, proteinCollection.proteome_id)
    longest_isoforms = []
    non_longest_isoforms = []
    all_isoforms = []
    for gene_id in data:
        longest_parsed = 0
        for protein_id, length in sorted(data[gene_id].items(), key=operator.itemgetter(1), reverse=True):
            if longest_parsed == 0:
                longest_isoforms.append(protein_id)
                longest_parsed = 1
            else:
                non_longest_isoforms.append(protein_id)
        all_isoforms.append("%s\t%s" % (gene_id, ",".join(str(protein_id) for (protein_id, length) in sorted(data[gene_id].items(), key=operator.itemgetter(1), reverse=True))))
    longest_isoforms = set(longest_isoforms)
    non_longest_isoforms = set(non_longest_isoforms)
    print "[+] %s proteinIDs in FASTA file" % (proteinCollection.protein_count)
    count_longest_isoforms = len(longest_isoforms)
    count_non_longest_isoforms = len(non_longest_isoforms)
    count_allocated = count_longest_isoforms + count_non_longest_isoforms
    fraction_fas_protein_ids_longest_isoform = count_longest_isoforms/proteinCollection.protein_count
    fraction_fas_protein_ids_non_longest_isoform = count_non_longest_isoforms/proteinCollection.protein_count
    fraction_fas_protein_ids_allocated = count_allocated/proteinCollection.protein_count

    if all_isoforms:
        print "[+] Writing %s all isoforms to %s ..." % (len(all_isoforms), out_all_isoforms_by_gene_f)
        with open(out_all_isoforms_by_gene_f, 'w') as out_all_isoforms_by_gene_fh:
            out_all_isoforms_by_gene_fh.write("\n".join(all_isoforms) + "\n")
    if longest_isoforms:
        print "[+] Writing %s (%s) proteinIDs to %s ..." % (count_longest_isoforms, "{:.1%}".format(fraction_fas_protein_ids_longest_isoform), out_longest_isoforms_f)
        with open(out_longest_isoforms_f, 'w') as out_longest_isoforms_fh:
            out_longest_isoforms_fh.write("\n".join(sorted(longest_isoforms)) + "\n")
    if non_longest_isoforms:
        print "[+] Writing %s (%s) proteinIDs to %s ..." % (count_non_longest_isoforms, "{:.1%}".format(fraction_fas_protein_ids_non_longest_isoform), out_non_longest_isoforms_f)
        with open(out_non_longest_isoforms_f, 'w') as out_non_longest_isoforms_fh:
            out_non_longest_isoforms_fh.write("\n".join(sorted(non_longest_isoforms)) + "\n")

    if fraction_fas_protein_ids_allocated == 1.0:
        if fraction_fas_protein_ids_longest_isoform == 1.0:
            print "[+] %s of proteins from FASTA allocated to %s" % ("{:.1%}".format(fraction_fas_protein_ids_allocated), out_longest_isoforms_f)
        else:
            print "[+] %s of proteins from FASTA allocated to %s and %s" % ("{:.1%}".format(fraction_fas_protein_ids_allocated), out_longest_isoforms_f, out_non_longest_isoforms_f)
        print "[+] Writing stats to %s ..." % (out_stats_f)
        with open(out_stats_f, 'w') as out_stats_fh:
            out_stats_fh.write("file=%s,total=%s,longest_isoforms=%s,longest_isoforms_perc=%s,non_longest_isoforms=%s,non_longest_isoforms_perc=%s,allocated_perc=%s,result=done\n" % (
                    fasta_f,
                    proteinCollection.protein_count,
                    count_longest_isoforms,
                    "{:.1%}".format(fraction_fas_protein_ids_longest_isoform),
                    count_non_longest_isoforms,
                    "{:.1%}".format(fraction_fas_protein_ids_non_longest_isoform),
                    "{:.1%}".format(fraction_fas_protein_ids_allocated)))




def parse_fasta(fasta_f, mapping_f, fs, fe, fdel):
    proteome_id = basename(fasta_f).split(".")[0]
    proteinCollection = ProteinCollection(proteome_id)
    i = 0
    print "## Debug ..."
    mapping_dict = {}
    if mapping_f:
        mapping_dict = parse_mapping_dict(mapping_f)
    for protein_id, protein_length in read_fasta_length(fasta_f):
        fas_protein_id = protein_id # this is used for output
        gff_protein_id = ''
        if mapping_dict:
            if fas_protein_id in mapping_dict:
                gff_protein_id = mapping_dict[fas_protein_id]
            else:
                print "[-] protein-ID %s not in mapping file" % (fas_protein_id)
        elif not fe:
            gff_protein_id = eval(fdel).join(protein_id.split(eval(fdel))[fs:]) # this is used for gff
        else:
            gff_protein_id = eval(fdel).join(protein_id.split(eval(fdel))[fs:fe]) # this is used for gff
        if i <= DEBUG_ENTRIES:
            print "# FASTA-ID=%s, GFF-ID=%s, LEN=%s" % (fas_protein_id, gff_protein_id, protein_length)
        elif i == DEBUG_ENTRIES + 1:
            print "# ..."
        else:
            pass
        i += 1
        proteinCollection.add_protein(fas_protein_id, gff_protein_id, protein_length)

    return proteinCollection

def parse_gff(gff_f, gff_type):

    protein_id_in_CDS_based_gff_types = set(["NCBI", "NHGRI", "SchistoDB", "OIST", "HGC", "OIST_MGU", "JGI", "VectorBase", "BeeBase", "Gnomon", "HGC", "ensembl_havana", "ensembl", "Ensembl_Metazoa", "BeetleBase", "NasoniaBase", "AphidBase", "ORCAE", "FlyBase"])
    gene_based_gff_types = set(["AUGUSTUS", "WormBase", "WormBase_imported", "transdecoder"])

    annotationCollection = AnnotationCollection()
    gene_id_by_protein_id = {}
    gene_id_by_mRNA_id = {}

    gff_protein_ids_with_mRNA = []
    gff_protein_ids_not_in_fas = []
    gff_protein_ids_in_fas = []

    if gff_type in protein_id_in_CDS_based_gff_types:
        for record_type, ninth_by_key in read_gff(gff_f):
            if record_type == 'CDS':
                mRNA_id = ninth_by_key['Parent'].replace("transcript:", "").replace("Transcript:", "")
                gff_protein_id = ninth_by_key.get('protein_id', None)
                if not gff_protein_id in gene_id_by_protein_id:
                    if mRNA_id in gene_id_by_mRNA_id:  # mRNA has been seen before/protein is real
                        gene_id = gene_id_by_mRNA_id[mRNA_id]
                        gene_id_by_protein_id[gff_protein_id] = gene_id
                        rescue_gff_protein_id = eval(fdel).join(gff_protein_id.split(eval(fdel))[0:-1])
                        if gff_protein_id in proteinCollection.fas_protein_id_by_gff_protein_id:
                            annotationCollection.add_protein(gene_id, gff_protein_id)
                            gff_protein_ids_in_fas.append(gff_protein_id)
                        elif rescue_gff_protein_id in proteinCollection.fas_protein_id_by_gff_protein_id:
                            print "# rescuing GFF-ID %s as %s" % (gff_protein_id, rescue_gff_protein_id)
                            annotationCollection.add_protein(gene_id, rescue_gff_protein_id)
                            gff_protein_ids_in_fas.append(rescue_gff_protein_id)
                        else:
                            gff_protein_ids_not_in_fas.append(gff_protein_id)
                        gff_protein_ids_with_mRNA.append(gff_protein_id)
            elif record_type == 'mRNA':
                gene_id = ninth_by_key['Parent'].replace("gene:", "").replace("Gene:", "")
                mRNA_id = ninth_by_key['ID'].replace("transcript:", "").replace("Transcript:", "")
                if gene_id in gene_id_by_mRNA_id:
                    sys.exit("[-] Collision of mRNA_id %s." % (mRNA_id))
                gene_id_by_mRNA_id[mRNA_id] = gene_id
            else:
                pass
    elif gff_type in gene_based_gff_types:
        for record_type, ninth_by_key in read_gff(gff_f):
            if record_type == 'mRNA' or record_type == "mrna": # don't ask ...
                gene_id = ninth_by_key['Parent'].replace("gene:", "").replace("Gene:", "")
                gff_protein_id = ninth_by_key['ID'].replace("transcript:", "").replace("Transcript:", "")
                rescue_gff_protein_id = eval(fdel).join(gff_protein_id.split(eval(fdel))[0:-1])
                if gff_protein_id in proteinCollection.fas_protein_id_by_gff_protein_id:
                    annotationCollection.add_protein(gene_id, gff_protein_id)
                    gff_protein_ids_in_fas.append(gff_protein_id)
                elif rescue_gff_protein_id in proteinCollection.fas_protein_id_by_gff_protein_id:
                    print "# rescuing GFF-ID %s as %s" % (gff_protein_id, rescue_gff_protein_id)
                    annotationCollection.add_protein(gene_id, rescue_gff_protein_id)
                    gff_protein_ids_in_fas.append(rescue_gff_protein_id)
                else:
                    gff_protein_ids_not_in_fas.append(gff_protein_id)
                gff_protein_ids_with_mRNA.append(gff_protein_id)
            else:
                pass
    else:
        sys.exit("[-] Unknown gff_type %s." % (gff_type))

    annotationCollection.gff_protein_ids_with_mRNA = set(gff_protein_ids_with_mRNA)
    annotationCollection.gff_protein_ids_in_fas = set(gff_protein_ids_in_fas)
    if len(annotationCollection.gff_protein_ids_in_fas) == 0:
        sys.exit("[-] None of the GFF protein IDs seems to be in FASTA\n\tExample of IDs in GFF:\n\t%s" % ",".join(gff_protein_ids_with_mRNA[0:10]))
    annotationCollection.gff_protein_ids_not_in_fas = set(gff_protein_ids_not_in_fas)
    return annotationCollection

class ProteinCollection():
    def __init__(self, proteome_id):
        self.proteome_id = proteome_id
        self.fas_protein_id_by_gff_protein_id = {}
        self.protein_length_by_gff_protein_id = {}
        self.protein_count = 0

    def add_protein(self, fas_protein_id, gff_protein_id, protein_length):
        self.fas_protein_id_by_gff_protein_id[gff_protein_id] = fas_protein_id
        self.protein_length_by_gff_protein_id[gff_protein_id] = protein_length
        self.protein_count += 1

class AnnotationCollection():
    def __init__(self):
        self.gff_protein_ids_by_gene_id = {}
        self.fas_protein_id_by_gff_protein_id = {}
        self.protein_length_by_gff_protein_id = {}

        self.gene_count = 0
        self.protein_count = 0
        self.gff_protein_ids_with_mRNA = set()
        self.gff_protein_ids_in_fas = set()
        self.gff_protein_ids_not_in_fas = set()
        self.fas_protein_ids = set()
        self.fas_protein_ids_in_gff_with_mRNA = set()
        self.fas_protein_ids_not_in_gff = set()


    def update(self, proteinCollection):
        for gff_protein_id, fas_protein_id in proteinCollection.fas_protein_id_by_gff_protein_id.items():
            if gff_protein_id in self.gff_protein_ids_with_mRNA:
                self.fas_protein_ids_in_gff_with_mRNA.add(fas_protein_id)
            else:
                self.fas_protein_ids_not_in_gff.add(fas_protein_id)
            self.fas_protein_ids.add(fas_protein_id)
        self.fas_protein_id_by_gff_protein_id = proteinCollection.fas_protein_id_by_gff_protein_id
        self.protein_length_by_gff_protein_id = proteinCollection.protein_length_by_gff_protein_id

    def write_stats(self):
        print "# gff_protein_ids_with_mRNA = %s" % (len(self.gff_protein_ids_with_mRNA))
        print "# gff_protein_ids_in_fas = %s" % (len(self.gff_protein_ids_in_fas))
        print "# gff_protein_ids_not_in_fas = %s" % (len(self.gff_protein_ids_not_in_fas))
        print "# fas_protein_ids = %s" % (len(self.fas_protein_ids))
        print "# fas_protein_ids_in_gff_with_mRNA = %s" % (len(self.fas_protein_ids_in_gff_with_mRNA))
        print "# fas_protein_ids_not_in_gff = %s" % (len(self.fas_protein_ids_not_in_gff))

    def get_longest_gff_protein_id(self, gene_id):
        lengths = []
        protein_ids = []
        for protein_id in self.gff_protein_ids_by_gene_id[gene_id]:
            if protein_id in proteinCollection.protein_length_by_gff_protein_id:
                protein_ids.append(protein_id)
                lengths.append(proteinCollection.protein_length_by_gff_protein_id[protein_id])
            else:
                sys.exit("[-] Unknown Protein ID %s" % (protein_id))
        if lengths:
            protein_idx, protein_length = max(enumerate(lengths), key=operator.itemgetter(1))
            return protein_ids[protein_idx], protein_length, list(set(protein_ids))
        else:
            return None, None, None

    def add_protein(self, gene_id, protein_id):
        if not gene_id in self.gff_protein_ids_by_gene_id:
            self.gff_protein_ids_by_gene_id[gene_id] = set()
            self.gene_count += 1
        self.gff_protein_ids_by_gene_id[gene_id].add(protein_id)
        self.protein_count += 1

def generate_output(out_prefix):
    out_non_longest_isoforms_f = "%s.non_longest_isoforms.txt" % (proteinCollection.proteome_id)
    out_longest_isoforms_f = "%s.longest_isoforms.txt" % (proteinCollection.proteome_id)
    out_fas_protein_ids_not_in_gff_f = "%s.fas_protein_ids_not_in_gff.txt" % (proteinCollection.proteome_id)
    out_all_isoforms_by_gene_f = "%s.isoforms_by_gene_id.txt" % (proteinCollection.proteome_id)
    out_stats_f = "%s.stats.txt" % (proteinCollection.proteome_id)
    if out_prefix:
        out_non_longest_isoforms_f = "%s.%s.non_longest_isoforms.txt" % (out_prefix, proteinCollection.proteome_id)
        out_longest_isoforms_f = "%s.%s.longest_isoforms.txt" % (out_prefix, proteinCollection.proteome_id)
        out_fas_protein_ids_not_in_gff_f = "%s.%s.out_fas_protein_ids_not_in_gff.txt" % (out_prefix, proteinCollection.proteome_id)
        out_all_isoforms_by_gene_f = "%s.%s.isoforms_by_gene_id.txt" % (out_prefix, proteinCollection.proteome_id)
        out_stats_f = "%s.%s.stats.txt" % (out_prefix, proteinCollection.proteome_id)

    longest_isoforms = []
    non_longest_isoforms = []
    all_isoforms = []
    for gene_id, gff_protein_ids in annotationCollection.gff_protein_ids_by_gene_id.items():
        if len(gff_protein_ids):
            longest_gff_protein_id, longest_protein_length, gff_protein_ids = annotationCollection.get_longest_gff_protein_id(gene_id)
            if longest_protein_length:
                all_isoforms_line = []
                longest_isoforms.append("%s" % (annotationCollection.fas_protein_id_by_gff_protein_id[longest_gff_protein_id]))
                for protein_id in gff_protein_ids:
                    all_isoforms_line.append(protein_id)
                    if not protein_id == longest_gff_protein_id:
                        non_longest_isoforms.append("%s" % (annotationCollection.fas_protein_id_by_gff_protein_id[protein_id]))
        all_isoforms.append("%s\t%s" % (gene_id, ",".join([annotationCollection.fas_protein_id_by_gff_protein_id[protein_id] for protein_id in gff_protein_ids])))
    longest_isoforms = set(longest_isoforms)
    non_longest_isoforms = set(non_longest_isoforms)
    print "[+] %s proteinIDs in FASTA file" % (proteinCollection.protein_count)
    count_longest_isoforms = len(longest_isoforms)
    count_non_longest_isoforms = len(non_longest_isoforms)
    count_allocated = count_longest_isoforms + count_non_longest_isoforms
    count_not_allocated = len(annotationCollection.fas_protein_ids_not_in_gff)
    fraction_fas_protein_ids_longest_isoform = count_longest_isoforms / proteinCollection.protein_count
    fraction_fas_protein_ids_non_longest_isoform = count_non_longest_isoforms / proteinCollection.protein_count
    fraction_fas_protein_ids_allocated = count_allocated / proteinCollection.protein_count
    fraction_fas_protein_ids_not_allocated = count_not_allocated / proteinCollection.protein_count

    if all_isoforms:
        print "[+] Writing %s all isoforms to %s ..." % (len(all_isoforms), out_all_isoforms_by_gene_f)
        with open(out_all_isoforms_by_gene_f, 'w') as out_all_isoforms_by_gene_fh:
            out_all_isoforms_by_gene_fh.write("\n".join(all_isoforms) + "\n")
    if longest_isoforms:
        print "[+] Writing %s (%s) proteinIDs to %s ..." % (count_longest_isoforms, "{:.1%}".format(fraction_fas_protein_ids_longest_isoform), out_longest_isoforms_f)
        with open(out_longest_isoforms_f, 'w') as out_longest_isoforms_fh:
            out_longest_isoforms_fh.write("\n".join(sorted(longest_isoforms)) + "\n")
    if non_longest_isoforms:
        print "[+] Writing %s (%s) proteinIDs to %s ..." % (count_non_longest_isoforms, "{:.1%}".format(fraction_fas_protein_ids_non_longest_isoform), out_non_longest_isoforms_f)
        with open(out_non_longest_isoforms_f, 'w') as out_non_longest_isoforms_fh:
            out_non_longest_isoforms_fh.write("\n".join(sorted(non_longest_isoforms)) + "\n")

    if fraction_fas_protein_ids_allocated == 1.0:
        if fraction_fas_protein_ids_longest_isoform == 1.0:
            print "[+] %s of proteins from FASTA allocated to %s" % ("{:.1%}".format(fraction_fas_protein_ids_allocated), out_longest_isoforms_f)
        else:
            print "[+] %s of proteins from FASTA allocated to %s and %s" % ("{:.1%}".format(fraction_fas_protein_ids_allocated), out_longest_isoforms_f, out_non_longest_isoforms_f)
        print "[+] Writing stats to %s ..." % (out_stats_f)
        with open(out_stats_f, 'w') as out_stats_fh:
            out_stats_fh.write("file=%s,total=%s,longest_isoforms=%s,longest_isoforms_perc=%s,non_longest_isoforms=%s,non_longest_isoforms_perc=%s,allocated_perc=%s,result=done\n" % (
                    fasta_f,
                    proteinCollection.protein_count,
                    count_longest_isoforms,
                    "{:.1%}".format(fraction_fas_protein_ids_longest_isoform),
                    count_non_longest_isoforms,
                    "{:.1%}".format(fraction_fas_protein_ids_non_longest_isoform),
                    "{:.1%}".format(fraction_fas_protein_ids_allocated)))
    else:
        print "[-] Only %s of proteins in FASTA allocated to %s and %s" % ("{:.1%}".format(fraction_fas_protein_ids_allocated), out_longest_isoforms_f, out_non_longest_isoforms_f)
        if annotationCollection.fas_protein_ids_not_in_gff:
            print "[+] Writing %s (%s) proteinIDs to %s ..." % (len(annotationCollection.fas_protein_ids_not_in_gff), "{:.1%}".format(fraction_fas_protein_ids_not_allocated), out_fas_protein_ids_not_in_gff_f)
            with open(out_fas_protein_ids_not_in_gff_f, 'w') as out_fas_protein_ids_not_in_gff_fh:
                out_fas_protein_ids_not_in_gff_fh.write("\n".join(annotationCollection.fas_protein_ids_not_in_gff) + "\n")
        print "[+] Writing stats to %s ..." % (out_stats_f)
        with open(out_stats_f, 'w') as out_stats_fh:
            out_stats_fh.write("file=%s,total=%s,longest_isoforms=%s,longest_isoforms_perc=%s,non_longest_isoforms=%s,non_longest_isoforms_perc=%s,allocated_perc=%s,result=check\n" % (
                    fasta_f,
                    proteinCollection.protein_count,
                    count_longest_isoforms,
                    "{:.1%}".format(fraction_fas_protein_ids_longest_isoform),
                    count_non_longest_isoforms,
                    "{:.1%}".format(fraction_fas_protein_ids_non_longest_isoform),
                    "{:.1%}".format(fraction_fas_protein_ids_allocated)))


if __name__ == "__main__":
    __version__ = 0.2
    args = docopt(__doc__)

    fasta_f = args['--fasta']
    fs = int(args['--fs'])
    try:
        fe = int(args['--fe'])
    except TypeError:
        fe = 0
    fdel = repr(args['--fdel'])
    print "# fs=%s, fe=%s, fdel=%s" % (fs, fe, fdel)
    gff_f = args['--gff3']
    mapping_f = args['--mapping']
    gff_type = args['--type']
    out_prefix = args['--outprefix']
    DEBUG_ENTRIES = 10
    proteinCollection = None
    print "[+] Start ..."

    if fasta_f and gff_f:
        print "[+] Parsing FASTA %s ..." % fasta_f
        proteinCollection = parse_fasta(fasta_f, mapping_f, fs, fe, fdel)

        print "[+] Parsing GFF3 %s ..." % gff_f
        annotationCollection = parse_gff(gff_f, gff_type)
        annotationCollection.update(proteinCollection)
        annotationCollection.write_stats()
        generate_output(out_prefix)
    elif fasta_f and gff_type == "custom":
        print "[+] Parsing FASTA %s ..." % fasta_f
        proteinCollection = parse_custom_fasta(fasta_f, fs, fe, fdel)
    else:
        print "[-] No input files given."
        sys.exit(__doc__.strip())

