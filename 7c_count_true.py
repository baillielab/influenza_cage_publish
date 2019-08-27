
infile = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr_gene_only_compared_all_true_times_two.tsv'

outfile = 'final_most_significant_all.tsv'
o = open(outfile, 'w')
o.write('gene\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\tpercent_snatched\tFDR\ttenmer\t2hv0h\t7hv2h\t24hv7h\t0h\t2h\t7h\t24\n')		


with open(infile, 'r') as inF:
	next(inF)
	for line in inF:
		counted = line.count('True')
		if counted >= 7:
			o.write(line)