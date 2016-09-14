#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from os.path import basename, isfile, abspath, splitext, join, exists
from os import getcwd, mkdir
from docopt import docopt, DocoptExit

def read_file(file_type, f):
    if not exists(f):
        kf_status.error(1)
    if file_type == "domain":


class DomainCollection():
    def __init__(self):
        self.f = None
        self.by_type_by_protein = None
        self.types = set()

    def parse_domains(self, domain_f):
        for line in read_file("domain", domain_f):

        self.f = domain_f


