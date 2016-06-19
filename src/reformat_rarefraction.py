#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import sys
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib as mat
mat.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("talk")
#sns.set(font='serif')
sns.set(style="whitegrid")

def reformat(infiles):
    for infile in infiles:
        inflation_value = ".".join(infile.split(".")[0:2]).lstrip("I")
        with open(infile) as fh:
            for line in fh:
                if not line.startswith("sample_size"):
                    temp = line.split()
                    print inflation_value + "\t" + temp[0] + "\t" + "\t".join(temp[4:])


if __name__ == "__main__":

    try:
        infiles = sys.argv[1:]
    except:
        sys.exit("")

    reformat(infiles)
