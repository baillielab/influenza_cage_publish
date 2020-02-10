from itertools import izip as zip, count
from collections import Counter
import os
import time
import re
import random
import logging
from functools import reduce
#from slacker import Slacker
import multiprocessing as mp
###############################################################################

import numpy as np 
from scipy.stats import stats

def pvalue_calculation(infile):

    sample = ((infile.split('_tenmers'))[0]).replace('donor', 'd')
    with open(allseqs_background) as inF:
        for line in inF:
            if 'seqs' in line:
                index_sample = ((line.strip()).split('\t')).index(sample)
            else:
                linea = line.split('\t')
                all_background = float(linea[index_sample])

    with open(allseqs_snatched) as inF:
        for line in inF:
            if 'seqs' in line:
                index_sample = ((line.strip()).split('\t')).index(sample)
            else:
                linea = line.split('\t')
                all_snatched = float(linea[index_sample])

    out = ('fishers_output_%s.txt')%(sample)
    o = open(out, 'w')
    with open(infile, 'r') as inF:
        for line in inF:
            linea = (line.strip()).split('\t')
            sequence = linea[0]
            snatch = float(linea[2])
            unsnatch = float(linea[1])

            A = snatch
            B = unsnatch
            C = all_snatched - snatch
            D = all_background - unsnatch

            oddsratio, pvalue = stats.fisher_exact([[A, B], [C,D]])

            outlist = [sequence, str(A), str(B), str(C), str(D), str(pvalue), str(oddsratio), '\n']
            output = '\t'.join(outlist)
            o.write(output)

###############################################################################

output_path = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
resource = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/'
tenmers = resource + '10mercount_tables/'

allseqs_background = tenmers + 'allseqs_background'
allseqs_snatched = tenmers + 'allseqs_snatched'

files = []
dirs = os.listdir(output_path)
for fle in dirs:
    if '.tsv' in fle:
        files.append(fle)

print files

pool = mp.Pool(processes=16)
results = pool.map(pvalue_calculation, files)

###############################################################################

#Fishers_function
'''
CONTINGENCY TABLE:
                snatched                background
this_seq        this_snatched   (A)     this_background     (B)
other_seqs      all - this      (C)     all - this          (D)
all_seqs        all_snatched            all_background

We use this table to find the p-value:

>>> import scipy.stats as stats
oddsratio, pvalue = stats.fisher_exact([[A, B], [C,D]])

'''
###############################################

        

###############################################