
yorn_file = 'timepoint_outfile_yorn.tsv'

yorn = {}

with open(yorn_file, 'r') as inF:
    for line in inF:
        linea = (line.strip()).split('\t')
        yorn[linea[0]] = linea[13] + '\t' + linea[14] + '\t' + linea[15]

print yorn['TCGGAGCGGT']



gene_dict = {}

name_dir = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/4_leaders/leader_tables/2_name/'

name_file = name_dir + 'promoter_list_filtered_iav_free.tsv'

with open(name_file, 'r') as inF:
    for line in inF:
        linea = (line.strip()).split('\t')
        if linea[-1] == 'NA':
            pass
        else:
            if linea[-1] in gene_dict:
                if linea[0] in gene_dict[linea[-1]]:
                    pass
                else:
                    gene_dict[linea[-1]].append(linea[0])
            else:
                gene_dict[linea[-1]] = [linea[0]]


print gene_dict['CTSB']

genes = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/4_leaders/leader_tables/4_genes/'

gene_list_file = genes + 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr_gene_only.tsv'
gene_list = {}

info_file = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr_gene_only_compared.tsv'
o = open(info_file, 'w')
o.write('gene\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\tpercent_snatched\tFDR\ttenmer\t2hv0h\t7hv2h\t24hv7h\n')
with open(gene_list_file, 'r') as inF:
    next(inF)
    for line in inF:
        line = line.strip()
        linea = line.split('\t')
        gene = linea[0]
        try:
            tenmer_list = gene_dict[gene]

            for tenmer in tenmer_list:
                output = yorn[tenmer]
                output = line + '\t' + tenmer + '\t' + output + '\n'
                #output = '\t'.join(output)
                o.write(output)
        except KeyError:
            output = 'NA\tNA\tNA\n'
            output = line + output
            o.write(output)


#compare

output = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr_gene_only_compared_all_true.tsv'
o = open(output, 'w')
o.write('gene\tsummed_or\tsummed_pvalue\tsnatch\tnotsnatch\tpercent_snatched\tFDR\ttenmer\t2hv0h\t7hv2h\t24hv7h\n')

with open(info_file, 'r') as inF:
    next(inF)
    for line in inF:
        counted = line.count('True')
        if counted >=3:
            o.write(line)


