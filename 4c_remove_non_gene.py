#remove 'chr'

infile = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr.tsv'
outfile = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr_gene_only.tsv'
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