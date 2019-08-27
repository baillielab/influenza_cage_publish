# import pandas
import pandas as pd
# import matplotlin
import matplotlib.pyplot as plt
# import seaborn
from collections import Counter
import seaborn as sns
import os
import math
import scipy.stats as ss
from statsmodels.stats.multicomp import pairwise_tukeyhsd
plt.style.use('ggplot')

#sns.set_style("ticks")

def Find_Files(path, string, exclude):
    import os
    files = []
    dirs = os.listdir(path)
    for fle in dirs:
        if string in fle:
            if exclude not in fle:
                if '.NA' not in fle:
                    files.append(fle)
    return files

def pvalue_recalc(x):
    if x == 0:
        return 300
    else:
        return float(-(math.log10(x)))

curr = os.getcwd()
files = Find_Files(curr, 'tsv.up.type', '.NA')
print files
for infile in files:

    out = infile + 'out'
    o = open(out, 'w')

    print infile
    types = []
    with open(infile, 'r') as inF:
        for line in inF:
            linea = line.split('\t')
            types.append(linea[-1])


    counted = Counter(types)
    o.write(str(counted))


    hour = (infile.split('_'))[3]

    df = pd.read_csv(infile, header = None, sep='\t')
    df.columns = ["sequence", "pvalue", "genename", "type"]
    df = df[df['type'] != 'antisense']
    df = df[df['type'] != 'processed_transcript']
    df = df[df['type'] != 'pseudo']
    df = df[df['type'] != 'polymorphic_pseudogene']
    df = df[df['type'] != 'sense_intronic']

    print df.head()



    # df = df[df['type'] != 'processed_transcript']
    df['pvalue'] = df['pvalue'].apply(lambda x: pvalue_recalc(x))

    #axes = sns.violinplot(x=df["type"],  y=df["pvalue"] )
    #plt.savefig(infile + '.pdf')

    sns.stripplot(x=df["type"],  y=df["pvalue"], palette="Blues", jitter=0.2, alpha = 0.5, edgecolor = 'grey', size= 5)
    sns.boxplot(x=df["type"],  y=df["pvalue"], palette="Blues", showfliers=False, boxprops=dict(alpha=0.7))#, jitter = 0.25)
    #sns.swarmplot(x=df["type"], y=df["pvalue"], color="grey")
    

    #sns.violinplot(x=df["type"],  y=df["pvalue"] )
    #sns.stripplot(x=df["type"],  y=df["pvalue"], jitter=0.25)

    #axes.get_figure().gca().set_title(hour)
    #axes.set_ylim(0, 350)
    #axes.set_axis_bgcolor('white')
    #axes.grid(False)

    p = plt.gca()
    p.set_xlabel("RNA type", fontsize =12)
    p.set_ylabel(("-log%s(FDR)")%('$_{10}$'), fontsize=12)
    plt.suptitle("")
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True) # labels along the bottom edge are off
    plt.title(('RNA type of snatched 10mers'))
    plt.savefig(infile + '.pdf')






    lincRNA = df[df['type'] == 'lincRNA']
    miRNA = df[df['type'] == 'miRNA']
    protein_coding = df[df['type'] == 'protein_coding']
    snRNA = df[df['type'] == 'snRNA']
    snoRNA = df[df['type'] == 'snoRNA']

    result = ss.f_oneway(lincRNA['pvalue'], miRNA['pvalue'], protein_coding['pvalue'], snRNA['pvalue'], snoRNA['pvalue'])
    print result


    from statsmodels.stats.multicomp  import pairwise_tukeyhsd
    from statsmodels.stats.multicomp  import MultiComparison
    import statsmodels.stats.multitest as multi

    mod = MultiComparison(df['pvalue'], df['type'])
    print mod.tukeyhsd()
    o.write(str(mod.tukeyhsd()))
    print mod.groupsunique
    o.write(str(mod.groupsunique))
