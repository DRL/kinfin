#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: nemnog_eval.py   -g <FILE> -n <FILE> -c <FILE>
                        [-h|--help]

    Options:
        -h --help                           show this
        -c, --category_file <FILE>          Category file
        -g, --groups <FILE>                 OrthologousGroups.txt produced by OrthoFinder
        -n, --nemnog <FILE>                 nemNOG members file

"""

from __future__ import division
from docopt import docopt
import sys
from os.path import basename, isfile, abspath, splitext, join, exists

class NemnogObj():
    def __init__(self, nemnog_cluster_id, nemnog_protein_count, nemnog_species_count, nemnog_functional_category, nemnog_proteins):
        self.cluster_id = nemnog_cluster_id
        self.protein_count = nemnog_protein_count
        self.protein_found = nemnog_protein_count
        self.species_count = nemnog_species_count
        self.functional_category = nemnog_functional_category
        self.proteins = nemnog_proteins
        self.orthogroups = None
        self.not_found = None
        self.not_found_fraction = None

    def get_orthogroups(self):
        orthogroups = {}
        not_found = []
        for protein in self.proteins:
            if protein in nemnogID_to_orthogroup:
                orthogroup = nemnogID_to_orthogroup[protein]
                if not orthogroup in orthogroups:
                    orthogroups[orthogroup] = []
                orthogroups[orthogroup].append(protein)
            else:
                not_found.append(protein)
                self.protein_found -= 1
        self.orthogroups = orthogroups
        self.not_found = not_found
        self.not_found_fraction = len(not_found)/self.protein_count

def parse_nemnog(nemnog_f):
    '''
    nemnog v4.5
    TaxonomicLevel|GroupName|ProteinCount|SpeciesCount|COGFunctionalCategory|ProteinIDs
    '''
    nemnogs = []
    with open(nemnog_f) as nemnog_fh:
        for line in nemnog_fh:
            temp = line.rstrip("\n").split()
            nemnog_cluster_id = temp[1]
            nemnog_protein_count = int(temp[2])
            nemnog_species_count = int(temp[3])
            nemnog_functional_category = temp[4]
            nemnog_proteins  = temp[5].split(",")
            nemnogs.append(NemnogObj(nemnog_cluster_id, nemnog_protein_count, nemnog_species_count, nemnog_functional_category, nemnog_proteins))
    return nemnogs

def parse_categories(category_f):
    prefix_to_taxid = {} # lists
    rankIDs = None
    with open(category_f) as category_fh:
        for l in category_fh:
            if l.startswith("#"):
                rankIDs = [x.strip() for x in l.lstrip("#").rstrip("\n").split(",")]
            else:
                if l.strip():
                    line = l.rstrip("\n").split(",")
                    if not len(line) == len(rankIDs):
                        sys.exit("[ERROR] - number of columns in line differs from header\n\t%s\n\t%s" % (rankIDs, line))
                    proteomeID = line[0]
                    for idx, levelID in enumerate(line):
                        rankID = rankIDs[idx]
                        if rankID == 'taxid':
                            if not levelID in prefix_to_taxid:
                                prefix_to_taxid[proteomeID] = []
                            prefix_to_taxid[proteomeID] = levelID
    return prefix_to_taxid

def parse_groups(groups_f):
    nemnogID_to_orthogroup = {}
    with open(groups_f) as group_fh:
        for line in group_fh:
            clusterID, protein_string = line.rstrip("\n").split(": ")
            for protein in protein_string.split():
                temp = protein.split(".")
                prefix_id = temp[0]
                translated_protein = None
                try:
                    taxid = prefix_to_taxid[prefix_id]
                    translated_protein = taxid + "." + ".".join(temp[1:])
                except:
                    translated_protein = protein
                nemnogID_to_orthogroup[translated_protein] = clusterID
    return nemnogID_to_orthogroup

def traverse_ids(nemnog_id, seen, path):
    nemnog_mapping_to_ogs = nemnog_fraction_of_orthogroup[nemnog_id].keys()
    for og_id in nemnog_mapping_to_ogs:
        path.append("%s=%s" % (og_id, nemnog_fraction_of_orthogroup[nemnog_id][og_id]))

        seen.add(og_id)
        nemnog_mapping_to_ogs = orthogroup_fraction_of_nemnog[og_id].keys()
        for nn_id in nemnog_mapping_to_ogs:
            if not nn_id in seen:
                path.append("%s=%s" % (nn_id, orthogroup_fraction_of_nemnog[og_id][nn_id]))
                seen.add(nn_id)
                seen, path = traverse_ids(nn_id, seen, path)
            else:
                pass
    return seen, path

def get_fractions():
    orthogroup_fraction_of_nemnog = {}
    nemnog_fraction_of_orthogroup = {}
    for nemnog in nemnogs:
        nemnog.get_orthogroups()
        for orthogroup, proteins in nemnog.orthogroups.items():
            if not orthogroup in orthogroup_fraction_of_nemnog:
                orthogroup_fraction_of_nemnog[orthogroup] = {}
            orthogroup_fraction_of_nemnog[orthogroup][nemnog.cluster_id] = len(proteins)/nemnog.protein_found
            if not nemnog.cluster_id in nemnog_fraction_of_orthogroup:
                nemnog_fraction_of_orthogroup[nemnog.cluster_id] = {}
            nemnog_fraction_of_orthogroup[nemnog.cluster_id][orthogroup] = len(proteins)/nemnog.protein_found
    return orthogroup_fraction_of_nemnog, nemnog_fraction_of_orthogroup

def get_counts():
    nemnogs_complete = {}
    nemnogs_complete_merged = {}
    nemnogs_split = {}
    nemnogs_split_merged = {}
    for cluster_id in orthogroup_fraction_of_nemnog:
        nemnog_id_count = len(orthogroup_fraction_of_nemnog[cluster_id].keys())
        for nemnog_id in orthogroup_fraction_of_nemnog[cluster_id]:
            nemnog_fraction = orthogroup_fraction_of_nemnog[cluster_id][nemnog_id]
            if nemnog_id_count == 1:
                if nemnog_fraction == 1.0:
                    nemnogs_complete[nemnog_id] = cluster_id
                else:
                    if not nemnog_id in nemnogs_split:
                        nemnogs_split[nemnog_id] = []
                    nemnogs_split[nemnog_id].append(cluster_id)
            else:
                if nemnog_fraction == 1.0:
                    nemnogs_complete_merged[nemnog_id] = cluster_id
                else:
                    if not nemnog_id in nemnogs_split_merged:
                        nemnogs_split_merged[nemnog_id] = []
                    nemnogs_split_merged[nemnog_id].append(cluster_id)
    return nemnogs_complete, nemnogs_complete_merged, nemnogs_split, nemnogs_split_merged

def write_stats(nemnogs_complete, nemnogs_complete_merged, nemnogs_split, nemnogs_split_merged):
    out_f = "%s.nemnog_eval.txt" % basename(groups_f)
    with open(out_f, 'w') as out_fh:
        for nemnog in nemnogs:
            if not nemnog.not_found_fraction == 1.0:
                if nemnog.cluster_id in nemnogs_complete:
                    out_fh.write("%s\tcomplete\t%s\t%s\t%s\t%s\n" % (nemnog.cluster_id, nemnog.functional_category, nemnog.species_count, nemnog.not_found_fraction, nemnogs_complete[nemnog.cluster_id]))
                elif nemnog.cluster_id in nemnogs_complete_merged:
                    out_fh.write("%s\tcomplete-merged\t%s\t%s\t%s\t%s\n" % (nemnog.cluster_id, nemnog.functional_category, nemnog.species_count, nemnog.not_found_fraction, nemnogs_complete_merged[nemnog.cluster_id]))
                elif nemnog.cluster_id in nemnogs_split_merged:
                    out_fh.write("%s\tsplit-merged\t%s\t%s\t%s\t%s\n" % (nemnog.cluster_id, nemnog.functional_category, nemnog.species_count, nemnog.not_found_fraction, ",".join(nemnogs_split_merged[nemnog.cluster_id])))
                elif nemnog.cluster_id in nemnogs_split:
                    out_fh.write("%s\tsplit\t%s\t%s\t%s\t%s\n" % (nemnog.cluster_id, nemnog.functional_category, nemnog.species_count, nemnog.not_found_fraction, ",".join(nemnogs_split[nemnog.cluster_id])))
                else:
                    sys.exit('[-] unknown id : %s' % nemnog.cluster_id)


    #with open(out_f, 'w') as out_fh:

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    try:
        groups_f = args['--groups']
        nemnog_f = args['--nemnog']
        category_f = args['--category_file']
    except docopt.DocoptExit:
        print __doc__.strip()

    print "[+] Starting analysis ..."
    print "[+] Parse categories %s ..." % category_f
    prefix_to_taxid = parse_categories(category_f)
    print "[+] Parse groups %s ..." % groups_f
    nemnogID_to_orthogroup = parse_groups(groups_f)
    print "[+] Parse nemNOGs %s ..." % nemnog_f
    nemnogs = parse_nemnog(nemnog_f)
    print "[+] Calculate nemnog-fractions of orthogroups ..."
    orthogroup_fraction_of_nemnog, nemnog_fraction_of_orthogroup = get_fractions()
    print "[+] Calculate counts ..."
    nemnogs_complete, nemnogs_complete_merged, nemnogs_split, nemnogs_split_merged = get_counts()
    print "[+] Writing files ..."
    write_stats(nemnogs_complete, nemnogs_complete_merged, nemnogs_split, nemnogs_split_merged)

'''
for file in ../data/OrthologousGroups_I*.txt ; do ~/git/clusterbuster/src/nemnog_eval.py -g $file -c ../data/SpeciesClassification.txt -n ../resources/nemNOG.members.tsv; done

By category : Totals
5751 S Function unknown (38.81%)
1641 T Signal transduction mechanisms (11.07%)
1008 O Posttranslational modification, protein turnover, chaperones
 862 U Intracellular trafficking, secretion, and vesicular transport
 826 K Transcription
 453 G Carbohydrate transport and metabolism
 440 P Inorganic ion transport and metabolism
 403 I Lipid transport and metabolism
 370 J Translation, ribosomal structure and biogenesis
 347 A RNA processing and modification
 299 Z Cytoskeleton
 299 C Energy production and conversion
 261 L Replication, recombination and repair
 253 W Extracellular structures
 223 E Amino acid transport and metabolism
'''
