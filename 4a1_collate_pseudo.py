
rnu = []
snrna = []

import os
import scipy.stats as stats
# from functools import reduce

total_background = 184630130
total_snatched = 8645162

infile = 'promoter_list_filtered_iav_free_summed_genes_5.tsv'
outfile = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo.tsv'
o = open(outfile, 'w')


with open(infile, 'r') as inF:
    for line in inF:
        if 'RNU' in line:
            linea = line.split('\t')
            if '-' in linea[0]:
                gene_name = (linea[0].split('-'))
                snrna.append(gene_name[0])
                rnu.append(line)
            elif 'RNU' in linea[0]:
                snrna.append(linea[0])
                rnu.append(line)
            else:
                o.write(line)
        else:
            o.write(line)

rnu_list =(list(set(snrna)))

rnu_list = (sorted(rnu_list, reverse=True))

i = 0
while i < len(rnu_list):
    gene = rnu_list[i]
    print gene
    snatched = 0
    not_snatched = 0

    for entry in rnu:
        if gene in entry:
            ind = rnu.index(entry)
            linea = entry.split('\t')
            snatched = snatched + float(linea[3])
            not_snatched = not_snatched + float(linea[4])
            del rnu[ind]

    a = snatched
    b = not_snatched
    c = total_snatched
    d = total_background

    oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
    output = [gene, str(oddsratio), str(pvalue), str(snatched), str(not_snatched), '\n']

    output = '\t'.join(output)
    o.write(output)

    i = i + 1


