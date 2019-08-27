#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

import os
import random
from itertools import izip as zip, count
from collections import Counter
import re
import multiprocessing as mp
import scipy.stats as stats
import statsmodels.stats.multitest as multi
import numpy as np
import pandas as pd
# from functools import reduce

def remove_influenza_from_tenmers(infile):
    influenza_6mer = ['GCAAAA', 'CAAAAG', 'AAAAGC', 'AAAGCA', 'AAGCAG', 'AGCAGG']

    outfile = (infile.split('.'))[0] + '_iav_free.tsv'
    o = open(outfile, 'w')

    line_list = []

    with open(infile, 'r') as inF:
        for line in inF:
            line_list.append(line)

    if len(set(line_list)) == len(line_list):
        for sequence in influenza_6mer:
            for line in line_list:
                if sequence in line:
                    linea = line.split('\t')
                    if sequence in linea[0]:
                        line_list.remove(line)

    print len(line_list)
    for line in line_list:
        o.write(line)

infile = 'promoter_list_filtered_iav_free.tsv'

remove_influenza_from_tenmers(infile)