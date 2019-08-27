infile = 'promoter_list_filtered_iav_free_summed_promoter_corr.tsv'

multi = 0

with open(infile, 'r') as inF:
	for line in inF:
		mm = line.count('@')
		if mm > 1:
			multi = multi + 1

print multi
