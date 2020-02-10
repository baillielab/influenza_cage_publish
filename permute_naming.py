#!/opt/local/bin/python
# -*- coding: UTF-8 -*-

#import os
import random
#from itertools import izip as zip, count
#from collections import Counter
#import re
import multiprocessing as mp
import scipy.stats as stats
import statsmodels.stats.multitest as multi
#import numpy as np
import pandas as pd
#import multiprocessing as mp
#from multiprocessing import Pool
# from functools import reduce
import sys

################################################

# FUNCTIONS AND CONSTANTS

################################################
def get_tenmers(infile):
    lst = []
    with open(infile, 'r') as inF:
        for line in inF:
            lst.append(line.strip())
    return lst

def ENTdictionary(sourcefile):
    ent_genename_dictionary = {}
    with open(sourcefile, 'r') as inF:
        next(inF)
        for line in inF:
            linea = line.split('\t')
            ent_genename_dictionary[linea[1]] = linea[2].strip()

    return ent_genename_dictionary

def remove_influenza_from_tenmers(infile):
    influenza_6mer = ['GCAAAA', 'CAAAAG', 'AAAAGC', 'AAAGCA', 'AAGCAG', 'AGCAGG']

    outfile = (infile.split('.'))[0] + '_iav_free.tsv'
    o = open(outfile, 'w')

    line_list = []

    with open(infile, 'r') as inF:
        for line in inF:
            line_list.append(line)

    if len(set(line_list)) == len(line_list):
        for line in line_list:
            i = 0
            y = 0
            while i < len(influenza_6mer):
                if influenza_6mer[i] in line:
                    y = y + 1
                    i = i + 1
                else:
                    i = i + 1
            if y == 0:
                o.write(line)

def fishers_for_promoter(input_file):
    values = (((input_file.split('.'))[0] + '_summed_promoter.tsv').split('/'))[-1]
    o = open(location_dir + values, 'w')
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

def correct_fishers(values):
    ids = []
    pvals = []
    perc = []

    with open(values, 'r') as inF:
        next(inF)
        for line in inF:
            linea = (line.strip()).split('\t')
            uniq_id = line.strip()
            pval = float(linea[2])
            ids.append(uniq_id)
            percent = (float(linea[3])/float(linea[4])) * 100

            perc.append(percent)
            pvals.append(pval)

    # r = rejected - bool
    # p = list of corrected pvals
    # s = sidak
    # b = bonf

    r, p, s, b = multi.multipletests(pvals,
                                     alpha=0.05,
                                     method='fdr_bh',
                                     is_sorted=False,
                                     returnsorted=False)
    p = list(p)

    values_out = (values.split('.'))[0] + '_corr.tsv'
    o = open(values_out, 'w')
    o.write('gene\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\tpercent_snatched\tFDR\n')

    o2 = open('pval.txt', 'w')
    o2.write('sidak pval = ' + str(s) + '\n')
    o2.write('bonf pval = ' + str(b) + '\n')
    o2.close()

    with open(values, 'r') as inF:
        next(inF)
        for line in inF:
            uniq_id = line.strip()
            if uniq_id in ids:
                ind = ids.index(uniq_id)
            else:
                print ('value not found %s') % (uniq_id)
                pass
            pvalue = p[ind]

            if float(pvalue) == 0.0:
                pvalue = 1e-350

            percent = perc[ind]
            line_out = line.strip() + '\t' + str(percent) + '\t' + str(pvalue) + '\n'
            o.write(line_out)

def TSSdictionary(infile):
    tss_dict = {}
    with open(infile, 'r') as inF:
        for line in inF:
            if '##' in line:
                pass
            else:
                line = line.split('\t')

                if ',' in line[1]:
                    prom = (line[1].split(','))[0]
                else:
                    prom = line[1]

                tss_dict[prom] = line[0]
    return tss_dict

def METAdictionary(tss_dict, infile):
    with open(infile, 'r') as inF:
        for line in inF:
            if '00' in line:
                pass
            else:
                line = line.split('\t')
                if ',' in line[1]:
                    prom = (line[1].split(','))[0]
                else:
                    prom = line[1]

                if prom not in tss_dict:
                    tss_dict[prom] = line[0]
                else:
                    pass
    return tss_dict

def assign_tss_to_promoter(infile):

    promoter_out = ((infile.split('.'))[0]) + '_ctss.tsv'
    o = open(promoter_out, 'w')
    o.write('gene\tctss\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\tpercent_snatched\tFDR\n')

    with open(infile, 'r') as inF:
        next(inF)
        for line in inF:
            linea = line.split('\t')
            prom = linea[0]
            if ',' in prom:
                prom = (prom.split(','))[0]
            if 'p@c' in prom:
                tss_out = prom
            else:
                try:
                    tss_out = tss_dict[prom]
                except KeyError:
                    tss_out = 'Not Found'
                    print prom + ' was not found in tss list? Checking all proms'

                    with open(infile, 'r') as inF:
                        for line in inF:
                            if '##' in line:
                                pass
                            else:
                                if prom in line:
                                    tss_out = (line.split('\t'))[0]
            outlist = [prom, tss_out] + linea[1:]
            output = '\t'.join(outlist)
            o.write(output)

def fishers_for_gene(input_file):
    values = (((input_file.split('.'))[0] + '_summed_genes.tsv').split('/'))[-1]
    o = open(values, 'w')
    o.write('gene\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\n')

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

    for gene in gene_set:
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
        o.write(outlist)

def collate_pseudo(infile):
    outfile = (infile.split('.'))[0] + '_pseudo.tsv'
    o = open(outfile, 'w')
    snrna = []
    rnu = []
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

    x = 0
    while x < len(rnu_list):
        gene = rnu_list[x]
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

        x = x + 1

def correct_fishers(values):
    ids = []
    pvals = []
    perc = []

    with open(values, 'r') as inF:
        next(inF)
        for line in inF:
            linea = (line.strip()).split('\t')
            uniq_id = line.strip()
            pval = float(linea[2])
            ids.append(uniq_id)
            percent = (float(linea[3])/float(linea[4])) * 100

            perc.append(percent)
            pvals.append(pval)

    # r = rejected - bool
    # p = list of corrected pvals
    # s = sidak
    # b = bonf

    r, p, s, b = multi.multipletests(pvals,
                                     alpha=0.05,
                                     method='fdr_bh',
                                     is_sorted=False,
                                     returnsorted=False)
    p = list(p)

    values_out = (values.split('.'))[0] + '_corr.tsv'
    o = open(values_out, 'w')
    o.write('gene\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\tpercent_snatched\tFDR\n')

    o2 = open('pval.txt', 'w')
    o2.write('sidak pval = ' + str(s) + '\n')
    o2.write('bonf pval = ' + str(b) + '\n')
    o2.close()

    with open(values, 'r') as inF:
        next(inF)
        for line in inF:
            uniq_id = line.strip()
            if uniq_id in ids:
                ind = ids.index(uniq_id)
            else:
                print ('value not found %s') % (uniq_id)
                pass
            pvalue = p[ind]

            if float(pvalue) == 0.0:
                pvalue = 1e-350

            percent = perc[ind]
            line_out = line.strip() + '\t' + str(percent) + '\t' + str(pvalue) + '\n'
            o.write(line_out)

def remove_non_gene(infile):
    outfile = (next_file.split('.'))[0] + '_gene_only.tsv'
    o = open(outfile, 'w')

    with open(infile, 'r') as inF:
        for line in inF:
            if 'chr' in line:
                if line[:3] == 'chr':
                    pass
                else:
                    o.write(line)
            else:
                o.write(line)

def pandas_get_top_1000(infile):
    df = pd.read_csv(infile, sep = '\t')
    df = df.sort_values('summed_or', ascending=False)

    df = df[df['FDR'] <= 0.05]
    df = df[df['summed_or'] > 0.0]
    df_top = df.head(1000)
    
    df.to_csv('all_tenmers.csv', sep = '\t', header = True)
    df_top.to_csv('all_tenmers_top_1000_snatched.csv', sep = '\t', header = True)
    permut_snatched = df_top['gene'].tolist()
    df = df.sort_values(by = 'summed_or', ascending=True)
    df_bot = df.head(1000)
    df_bot.to_csv('all_tenmers_top_1000_unsnatched.csv', sep = '\t', header = True)
    permut_unsnatched = df_bot['gene'].tolist()
    return permut_snatched, permut_unsnatched
    
    #top 100 from each group
snatched_f = 'snatched.txt'
unsnatched_f = 'unsnatched.txt'
snatched = get_tenmers(snatched_f)
unsnatched = get_tenmers(unsnatched_f)
# print 'SNATCHED'
# print snatched[0:100]
#now have 2 lists of 1000
#these are constants for this dataset
total_background = 184630130
total_snatched = 8645162
#these are constants for this dataset
gene_names = ENTdictionary('grch37_180401_geneid_transid_genename_mart_export.txt')
tss_dict = TSSdictionary('hg19.cage_peak_phase1and2combined_ann.txt')
meta_dict = METAdictionary(tss_dict, 'METADATA_U22_ENHANCERS_reformat.tsv')

print meta_dict['p3@HPSE2']

################################################

# FUNCTIONS

################################################



infile = sys.argv[1]#'promoter_list_unfiltered.tsv'
infile_outnum = ((infile.split('_'))[3]).replace('.tsv', '')

#infile_outnum = 'total'
score_file = (('score_file_%s.tsv')%(infile_outnum))
s = open(score_file, 'ab+')
s.write('snatched_score\tunsnatched_score\n')
s.close()

'''
i = 0
while i <= 1:
    #create promoter list file
    check_list = []
    #remove_output
    outfile = (('promoter_list_filtered_%s.tsv')%(str(i)))
    o = open(outfile, 'w')
    outlist = ['tenmer', 'hour', 'donor', 
               'FDR', 'OR', 'snatched', 'unsnatched', 
               'promoter', 'single_prom', 'gene', '\n']
    output = '\t'.join(outlist)
    o.write(output)
    with open(infile, 'r') as inF:
        next(inF)
        for line in inF:
            linea = (line.strip()).split('\t')
            if len(linea) >= 8:
                if '..' in linea[7]:
                    outprom = linea[7]
                    outgene = (outprom.split('@'))[1]
                    outlist = linea[0:7] + [outprom, outgene, '\n']
                else:

                    promoter = linea[7].split(',')
                    if len(promoter) > 1:
                        single_prom = random.choice(promoter) #randomly choose a promoter
                    else:
                        single_prom = promoter[0]
                    gene = (single_prom.split('@'))[1]

                    if 'ENST' in gene:
                        if gene in gene_names:
                            gene = gene_names[gene]

                    outprom = (single_prom.split('@'))[0] + '@' + gene
                    outgene = gene
                    outlist = linea[0:7] + [outprom, outgene, '\n']
            
            else:
                outlist = linea[0:7] + ['NA', 'NA', '\n']
            
            output = '\t'.join(outlist)
            o.write(output)
    
    o.close()
    infile_flu_remove = outfile
    '''
remove_influenza_from_tenmers(infile) #needs to be indented once if above code in use
#next_file ='iav_free_' + infile
next_file = ((infile.split('.'))[0])+ '_iav_free.tsv'

# promoter
# fishers_for_promoter(next_file)
# next_file = (next_file.split('.'))[0] + '_summed_promoter.tsv'
# correct_fishers(next_file)
# next_file = (next_file.split('.'))[0] + '_corr.tsv'
# assign_tss_to_promoter(next_file)

# gene
fishers_for_gene(next_file)
next_file = ((next_file.split('.'))[0])+ '_summed_genes.tsv'

collate_pseudo(next_file)
next_file = ((next_file.split('.'))[0])+ '_pseudo.tsv'

correct_fishers(next_file)
next_file = ((next_file.split('.'))[0])+ '_corr.tsv'

remove_non_gene(next_file)
next_file = ((next_file.split('.'))[0])+ '_gene_only.tsv'

perm_snatch, perm_unsnatch = pandas_get_top_1000(next_file)

snatch_score = len(set(perm_snatch) & set(snatched))
unsnatch_score = len(set(perm_unsnatch) & set(unsnatched))

print infile
print snatch_score
print unsnatch_score

s = open(score_file, 'a+')
s.write(str(snatch_score) + '\t' + str(unsnatch_score) + '\n')
s.close()

    #i = i + 1
