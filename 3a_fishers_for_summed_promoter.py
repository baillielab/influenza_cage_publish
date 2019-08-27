#!/opt/local/bin/python
# -*- coding: UTF-8 -*-
import os
import scipy.stats as stats

total_background = 184630130
total_snatched = 8645162

# from functools import reduce
capsnatch = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
location_dir = capsnatch + '4_leaders/leader_tables/2_name/'
infile = 'promoter_list_filtered_iav_free.tsv'
promoters = capsnatch + '4_leaders/leader_tables/3_promoters/'

def fishers_for_promoter(input_file):
    values = (((input_file.split('.'))[0] + '_summed_promoter.tsv').split('/'))[-1]
    promoters = capsnatch + '4_leaders/leader_tables/3_promoters/'
    o = open(promoters + values, 'w')
    o.write('promoter\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\n')

    promoters = []
    with open(input_file, 'r') as inF:
        next(inF)
        for line in inF:
            linea = line.split('\t')
            promoter = linea[7]
            promoters.append(promoter)

    promoter_set = list(set(promoters))

    for promoter in promoter_set:
        snatched = 0
        not_snatched = 0
        with open(input_file, 'r') as inF:
            next(inF)
            for line in inF:
                if promoter in line:
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
        output = [promoter, str(oddsratio), str(pvalue),
                  str(snatched), str(not_snatched), '\n']
        outlist = '\t'.join(output)
        o.write(outlist)



fishers_for_promoter(location_dir + infile)
