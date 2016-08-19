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

if __name__ == "__main__":
    contingency = [
        [100, 10000],
        [1, 10000000]
    ]
    print contingency
    calculate_fishers_exact_test(contingency)
    chisquare(contingency)
