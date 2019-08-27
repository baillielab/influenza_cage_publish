import pandas as pd

infile = 'promoter_list_filtered_iav_free_summed_genes_5_pseudo_corr_gene_only.tsv'

df = pd.read_csv(infile, sep = '\t')
df = df.sort('summed_or', ascending=False)
df = df[df['FDR'] <= 0.05]
df = df[df['summed_or'] > 0.0]

df_top = df.head(1000)

df_top.to_csv('top_1000_genes.tsv', sep = '\t', index = False)


df = df.sort('summed_or', ascending=True)

df_bot = df.head(1000)

df_bot.to_csv('bottom_1000_genes.tsv', sep = '\t', index = False)

df = df.sort('percent_snatched', ascending=False)

df_perc = df.head(1000)

df_perc.to_csv('top_1000_genes_perc.tsv', sep = '\t', index = False)