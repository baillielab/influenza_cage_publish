import os
import scipy.stats as stats
import statsmodels.stats.multitest as multi



infile = 'promoter_list_filtered_iav_free_summed_promoter.tsv'

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

correct_fishers(infile)