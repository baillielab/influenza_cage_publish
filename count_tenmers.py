from itertools import izip as zip, count
from collections import Counter
import os
import time
import re
import logging
from functools import reduce
#from slacker import Slacker
import multiprocessing as mp
###############################################################################
# get locations
###############################################################################
output_path = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
path = '/mnt/ris-fas1a/linux_groups2/fantom5/Clustering/FLU_TIMECOURSE/bamfiles/'
leader_file = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/tenmer_SC.txt'

#leader_file = 'tenmer_testfile.txt'
###############################################################################
query_list = []
query_flu_list = []
file_list = []
with open(leader_file, 'r') as inF:
    next(inF)
    for line in inF:
        line_split = line.strip()
        tenmer = line_split
        query_list.append(tenmer)
        #query_flu_list.append(tenmer + 'GCAAAAGCAGG')

dirs = os.listdir(path)
for fle in dirs:
    if '.sam' in fle:
        if 'mock' not in fle:
            file_list.append(fle)

###############################################################################
# locate_leaders
###############################################################################
def locate_leaders(filename):
    #filename = path + filename
    bamfiles_list = []
    bamfiles_list2 = []

    print filename

    with open(filename, 'r') as inF:
        for line in inF:
            if 'VHE' in line:
                line = line.split('\t')
                bamfiles_list.append((line[9])[0:10])
                #bamfiles_list2.append((line[9])[0:21])



    donor = (re.split('%20|%3a|%29', filename))[8]
    sample = (re.split('%20|%3a|%29', filename))[10]
    output = (output_path + '/%s_%s_tenmers_SC.tsv') % (sample, donor)

    segments = output_path + '3_segment/publish/segment/'
    segmented = segments + '20170718_GCAAAGCAGG_segmented.txt'

    flu_search = ('GCAAAAGCAGG_' + '%s_%s') % (sample, donor)
    flu_search = flu_search.replace('onor', '')

    query_flu_list = []
    with open(segmented, 'r') as inF:
        next(inF)
        for line in inF:
            if flu_search in line:
                leader = (line.split('\t'))[3]
                if len(leader) < 10:
                    pass
                else:
                    if len(leader) >= 10:
                        query_flu_list.append(leader[:10])


    #print (('%s has %s entries') % (filename, len(bamfiles_list2)))

    f = open(output, 'r')
    lines = f.readlines()
    f.close()

    length = len(lines)

    o = open(output, 'ab+')
    for i in range(length, len(query_list)):
        total = bamfiles_list.count(query_list[i])
        total_flu = bamfiles_list.count(query_flu_list[i])
        o.write(str(query_list[i]) + '\t' + str(total) + '\t' + str(total_flu))
        o.write('\n')

  
###############################################################################
# this part take s along time and the time differs between samfile.

pool = mp.Pool(processes=16)
results = pool.map(locate_leaders, file_list)

#locate_leaders('Monocyte-derived%20macrophages%20response%20to%20udorn%20influenza%20infection%2c%2002hr00min%2c%20donor2%20%28150_120%3aUd_2h%29.CNhs13647.13318-143A6.hg19.nobarcode.sam')
###############################################################################
'''
###############################################
import numpy as np 
from scipy.stats import stats

this_seq_snatched = 'thisseq_snatched'
this_seq_notsnatched = 'thisseq_notsnatched'
allseqs_background = 'allseqs_background'
allseqs_snatched = 'allseqs_snatched'

resource = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/'
tenmers = resource + '10mercount_tables/'
###############################################
#create lists with tenmers and samples required
samples = []
tenmers = []
total_background = []
total_snatched = []

files = []
dirs = os.listdir(output_path)
for fle in dirs:
    if '.tsv' in fle:
        files.append(fle)

print files


with open(this_seq_snatched) as inF:
    for line in inF:
        linea = (line.strip()).split('\t')
        if 'seqs' in line:
            for sample in linea[1:]:
                samples.append(sample)
        else:
            tenmers.append(linea[0])

with open(allseqs_background) as inF:
    next(inF)
    for line in inF:
        linea = (line.strip()).split('\t')
        for count in linea[1:]:
            total_background.append(float(count))

with open(allseqs_snatched) as inF:
    next(inF)
    for line in inF:
        linea = (line.strip()).split('\t')
        for count in linea[1:]:
            total_snatched.append(float(count))

###############################################
#Fishers_function



CONTINGENCY TABLE:
                snatched                background
this_seq        this_snatched   (A)     this_background     (B)
other_seqs      all - this      (C)     all - this          (D)
all_seqs        all_snatched            all_background

We use this table to find the p-value:

>>> import scipy.stats as stats
oddsratio, pvalue = stats.fisher_exact([[A, B], [C,D]])



###############################################

def Get_Value(infile, tenmer):
    with open(infile, 'r') as inF:
            for line in inF:
                if tenmer in line:
                    line = line.split('\t')
                    return float(line[i+1])

#write output
out = 'fishers_output.txt'
o = open(out, 'w')
o.write('##sample\t')
for sample in samples:
    o.write(sample + '\t')

#fishers exact on tenmers
for tenmer in tenmers:
    o.write('\n' + tenmer)
    i = 0
    while i < len(samples):
        A = Get_Value(this_seq_snatched, tenmer)
        B = Get_Value(this_seq_notsnatched, tenmer)
        C = total_snatched[i] - A
        D = total_background[i] - B
#perform statistical test
        oddsratio, pvalue = stats.fisher_exact([[A, B], [C,D]])
        o.write('\t' + str(oddsratio) + '|' + str(pvalue))
        i = i + 1
        print oddsratio, pvalue

###############################################
'''
