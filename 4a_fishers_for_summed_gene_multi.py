#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
import os
import scipy.stats as stats
import multiprocessing as mp
from multiprocessing import Pool

# from functools import reduce

total_background = 184630130
total_snatched = 8645162

def fishers_for_gene(gene):
    snatched = 0
    not_snatched = 0
    with open(input_file, 'r') as inF:
        next(inF)
        for line in inF:
            if gene in line:
                linea = (line.strip()).split('\t')
                this_snatched = float(linea[5])
                snatched = snatched + this_snatched
                this_not_snatched = float(linea[6])
                not_snatched = not_snatched + this_not_snatched
    a = snatched
    b = not_snatched
    c = total_snatched
    d = total_background

    oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
    output = [gene, str(oddsratio), str(pvalue),
              str(snatched), str(not_snatched), '\n']
    outlist = '\t'.join(output)
    o = open(values, 'ab+')
    o.write(str(outlist))
    o.close()



#def fishers_for_gene_wrap(input_file):

capsnatch = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
location_dir = capsnatch + '4_leaders/leader_tables/2_name/'
infile = 'promoter_list_filtered_iav_free.tsv'
genes = capsnatch + '4_leaders/leader_tables/4_genes/'
input_file = location_dir + infile
values = (((input_file.split('.'))[0] + '_summed_genes_5.tsv').split('/'))[-1]
#genes = capsnatch + '4_leaders/leader_tables/4_genes/'
o = open(genes + values, 'w')
o.write('gene\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\n')
o.close()
genes = []
with open(input_file, 'r') as inF:
    next(inF)
    for line in inF:
        linea = line.split('\t')
        gene = linea[8]
        if ',p' in gene:
            gene = (gene.split(',p'))[0]
        genes.append(gene)

gene_set = list(set(genes))

print len(gene_set)
pool = mp.Pool(processes=10)
pool.map(fishers_for_gene, gene_set)

#fishers_for_gene_wrap(location_dir + infile)
