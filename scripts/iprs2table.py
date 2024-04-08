#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage: ips_to_table.py      -i <FILE> [-o <STR>]
                            [--domain_sources <STRING>]
                            [--only_domain_go]
                            [--only_domain_ipr]
                            [-h|--help]

    Options:
        -h --help                               show this
        -i, --interproscan_f <FILE>             Interproscan file
        -o, --out_prefix <STR>                  Outprefix (default: )
        --domain_sources <STRING>               Collect information for those domain_sources [default: SignalP_EUK,Pfam]
                                                    - GO-terms and IPR-IDs domain counts are inferred for these
"""
import sys
import collections
from os.path import exists

import_errors = []
try:
    from docopt import docopt
except ImportError:
    import_errors.append("[ERROR] : Module \'Docopt\' was not found. Please install \'Docopt\' using \'pip install docopt\'")
if import_errors:
    sys.exit("\n".join(import_errors))


def read_file(f):
    if not f or not exists(f):
        sys.exit("[ERROR] - File '%s' does not exist." % (f))
    with open(f) as fh:
        for line in fh:
            yield line.rstrip("\n")


def parse_interproscan(inputObj):
    domain_sources = inputObj.domain_sources
    domains_ids_by_domain_source_by_protein_id = {}
    for line in read_file(inputObj.interproscan_f):
        temp = line.split()
        protein_id = temp[0]
        domain_source = temp[3]
        domain_id = temp[4]
        putative_ipr = [field for field in temp[6:] if field.startswith("IPR")]
        putative_gos = temp[-1]
        if protein_id not in domains_ids_by_domain_source_by_protein_id:
            domains_ids_by_domain_source_by_protein_id[protein_id] = {}
        if domain_source in domain_sources:
            if domain_source not in domains_ids_by_domain_source_by_protein_id[protein_id]:
                domains_ids_by_domain_source_by_protein_id[protein_id][domain_source] = []
            domains_ids_by_domain_source_by_protein_id[protein_id][domain_source].append(domain_id)
            if putative_gos.startswith("GO:"):
                go_terms = putative_gos.split("|")
                for go_term in go_terms:
                    if 'GO' not in domains_ids_by_domain_source_by_protein_id[protein_id]:
                        domains_ids_by_domain_source_by_protein_id[protein_id]['GO'] = set()
                    domains_ids_by_domain_source_by_protein_id[protein_id]['GO'].add(go_term)
            for ipr in putative_ipr:
                if "IPR" not in domains_ids_by_domain_source_by_protein_id[protein_id]:
                    domains_ids_by_domain_source_by_protein_id[protein_id]['IPR'] = []
                domains_ids_by_domain_source_by_protein_id[protein_id]['IPR'].append(ipr)
    return domains_ids_by_domain_source_by_protein_id


def output(inputObj, domains_ids_by_domain_source_by_protein_id):
    domain_sources = inputObj.domain_sources
    out_prefix = inputObj.out_prefix
    out_f = ''
    if out_prefix:
        out_f = "%s.functional_annotation.txt" % (out_prefix)
    else:
        out_f = "functional_annotation.txt"
    print("[+] Writing output to %s ..." % (out_f))
    output = []
    header = ["#protein_id"] + domain_sources
    output.append("\t".join(header))
    protein_id_counts_by_domain_source = {}
    unique_domain_ids_by_domain_source = {}
    for protein_id in sorted(domains_ids_by_domain_source_by_protein_id):
        line = []
        line.append(protein_id)
        for domain_source in domain_sources:
            if domain_source in domains_ids_by_domain_source_by_protein_id[protein_id]:
                if domain_source == "GO":
                    line.append(";".join(sorted(list(domains_ids_by_domain_source_by_protein_id[protein_id][domain_source]))))
                else:
                    domain_counts = collections.Counter(domains_ids_by_domain_source_by_protein_id[protein_id][domain_source])
                    # print ["%s:%s" % (domain_id, count) for domain_id, count in domain_counts.items()]
                    line.append(";".join(["%s:%s" % (domain_id, count) for domain_id, count in domain_counts.items()]))
                protein_id_counts_by_domain_source[domain_source] = protein_id_counts_by_domain_source.get(domain_source, 0) + 1
                if domain_source not in unique_domain_ids_by_domain_source:
                    unique_domain_ids_by_domain_source[domain_source] = set()
                for domain_id in domains_ids_by_domain_source_by_protein_id[protein_id][domain_source]:
                    unique_domain_ids_by_domain_source[domain_source].add(domain_id)
            else:
                line.append("None")
        output.append("\t".join(line))
    with open(out_f, "w") as out_fh:
        out_fh.write("\n".join(output) + "\n")
    counts = []
    for domain in domain_sources:
        counts.append("[=]\t%s = %s unique IDs in %s proteins" % (domain, len(unique_domain_ids_by_domain_source.get(domain, {})), protein_id_counts_by_domain_source.get(domain, 0)))
    print("\n".join(counts))


class InputObj():
    def __init__(self, args):
        # print args
        self.interproscan_f = args['--interproscan_f']
        self.out_prefix = args['--out_prefix']
        self.domain_sources = ["GO", "IPR"] + args['--domain_sources'].split(",")


if __name__ == "__main__":
    args = docopt(__doc__)
    print("[+] Start... ")
    inputObj = InputObj(args)
    print("[+] Parsing domain-IDs for %s" % (",".join(inputObj.domain_sources)))
    print("[+] Parsing domains in %s ..." % (inputObj.interproscan_f))
    domains_ids_by_domain_source_by_protein_id = parse_interproscan(inputObj)
    output(inputObj, domains_ids_by_domain_source_by_protein_id)

'''
Pending:
    - make sure that at least one domain is parsed (otherwise there is no GO-Terms)
'''
