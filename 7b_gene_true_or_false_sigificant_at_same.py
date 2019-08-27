
file_name = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr_gene_only_compared_all_true.tsv'


tenmer_list = []

with open(file_name, 'r') as inF:
	for line in inF:
		linea = line.split('\t')
		tenmer_list.append(linea[7])


times = ['0h', '2h', '7h', '24h']

file_input = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/4_leaders/leader_tables/1_tenmers/tenmer_table_full_corrected.tsv'

output = 'out.tsv'
o = open(output, 'w')


for tenmer in tenmer_list:
	hr0 = []
	hr2 = []
	hr7 = []
	hr24 = []

	hours = [hr0, hr2, hr7, hr24]

	with open(file_input, 'r') as inF:
		next(inF)
		for line in inF:
			if tenmer in line:
				linea = line.split('\t')
				ind = times.index(linea[1])
				hours[ind].append(linea[4])
	
	out = []

	for x in hours:
		if x.count('NA') >=3:
			out.append('False')
		else:
			out.append('True')
	outlist = '\t'.join([tenmer] + out + ['\n'])
	o.write(outlist)

o.close()

tenmer_dict = {}

with open(output, 'r') as inF:
	for line in inF:
		linea = line.split('\t')
		print linea[0]
		tenmer_dict[linea[0]] = linea[1] + '\t' + linea[2] + '\t' + linea[3] + '\t' + linea[4]

print tenmer_dict['ATACTTACCT']

out = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr_gene_only_compared_all_true_times_two.tsv'
o = open(out, 'w')
o.write('gene\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\tpercent_snatched\tFDR\ttenmer\t2hv0h\t7hv2h\t24hv7h\t0h\t2h\t7h\t24\n')		

with open(file_name, 'r') as inF:
	next(inF)
	for line in inF:
		line = line.strip()
		linea = line.split('\t')
		tenmer = linea[7]
		print tenmer
		try:
			is_true = tenmer_dict[tenmer]
			o.write(line +'\t' + is_true + '\n')
		except KeyError:
			pass
