#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from scipy import arange, stats
from scipy.stats import chi2_contingency
import numpy as np

def calculate_fishers_exact_test(list_of_lists):
    oddsratio, pvalue = stats.fisher_exact(list_of_lists)
    print "[+] Fisher's Exact Test"
    print oddsratio
    print pvalue

def chisquare(list_of_lists):
    obs = np.array(list_of_lists)
    g, p, dof, expctd = chi2_contingency(obs, lambda_="log-likelihood")
    print "[+] Chi-Square test (G-test)"
    print g
    print p

def parse_domains(domain_f):
    domains_by_proteinID = {}
    print "[STATUS] - Parsing domains from %s" % (domain_f)
    with open(domain_f) as fh:
        for line in fh:
            temp = line.rstrip("\n").split()
            proteinID = temp[0]
            domain = temp[1]
            domain_type = temp[3]
            domain_id = temp[4]
            evalue = temp[-3]
            stop = temp[-4]
            start = temp[-5]
            if len(" ".join(temp[5:-5])):
                desc = "\"%s\"" % " ".join(temp[5:-5])
            else:
                desc = None
            print proteinID, domain, domain_type, domain_id, evalue, start, stop, desc

if __name__ == "__main__":
    contingency = [
        [100, 10000],
        [1, 10000000]
    ]
    #print contingency
    #calculate_fishers_exact_test(contingency)
    #chisquare(contingency)

    parse_domains(sys.argv[1])
