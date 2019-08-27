import os
from collections import Counter

infile = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/resource/parsed_human_genes_all.txt'
f = open(infile, 'r')
lines = f.readlines()
f.close()

#get_type
get_type =  {}
i = 0
while i < len(lines):
    if '>' in lines[i]:
        name = ((lines[i].split('\t'))[4]).replace('>', '')
        type_ = ((lines[i].split('\t'))[-3])
        if '-' in name:
            name = (name.split('-'))[0] 
        get_type[name] = type_
    i = i + 1


def Find_Files(path, string):
    import os
    files = []
    dirs = os.listdir(path)
    for fle in dirs:
        if string in fle:
                files.append(fle)
    return files


def ENTdictionary(sourcefile):
    ent_genename_dictionary = {}
    with open(sourcefile, 'r') as inF:
        next(inF)
        for line in inF:
            linea = line.split('\t')
            ent_genename_dictionary[linea[1]] = linea[2]

    return ent_genename_dictionary

def pseudolist(sourcefile):
    pseudo_list = []
    with open(sourcefile, 'r') as inF:
        for line in inF:
            line = line.split('\t')
            pseudo_list.append(line[2])
    return pseudo_list

pseudo = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/resource/human_pseudogenes.txt'
pseudolist = pseudolist(pseudo)
gene_names = ENTdictionary('/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/resource/grch37_180401_geneid_transid_genename_mart_export.txt')
current = os.getcwd() + '/'
take_out = ['antisense', 'processed transcript']
bits = [",", "@", "-"]

fle = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/4_leaders/leader_tables/4_genes/promoter_list_filtered_iav_free_summed_genes_pseudo_corr_gene_only.tsv'

#files = ['summary_table_genes_adjusted.tsv']


output = str(fle) + '.up.type'
output3 = str(fle) + '.down.type'
output2 = str(fle) + '.type.NA'
o = open(output, 'w')
o2 = open(output2, 'w')
o3 = open(output3, 'w')

type_list = []
with open(fle, 'r') as inF:
    next(inF)
    for line in inF:
        if 'NaN' in line:
            pass
        else:
            linea = line.split('\t')
            seq = 'seq' #linea[0]
            name = linea[0]
            or_val = float(linea[1])
            pval = linea[6].strip()

            print or_val

            if or_val < 1:
                out = o3
            elif or_val == 0:
                pass
            else:
                print or_val
                out = o
            if ',' in name:
                name = (name.split(','))[0]
            if '-' in name:
                name = (name.split('-'))[0]
            if '@' in name:
                name = (name.split('@'))[1]
            if 'ENST' in name:
                try:
                    name = gene_names[name]
                except KeyError:
                    name = name
            name = name.strip()
            try:
                type_ = get_type[name]
                type_list.append(type_)
                for i in take_out:
                    if type_ == i:
                        outlist = [seq, pval, name, 'NA']#, loc]
                        output = '\t'.join(outlist)
                        o2.write(output + '\n')

                outlist = [seq, pval, name, type_,]# loc]
                output = '\t'.join(outlist)
                out.write(output + '\n')
            except KeyError:
                if 'RNU' in name:
                    outlist = [seq, pval, name, 'snRNA']#, loc]
                    output = '\t'.join(outlist)
                    out.write(output + '\n')
                elif 'RP11' in name:
                    outlist = [seq, pval, name, 'pseudo']#, loc]
                    output = '\t'.join(outlist)
                    out.write(output + '\n')
                elif 'RP13' in name:
                    outlist = [seq, pval, name, 'pseudo']#, loc]
                    output = '\t'.join(outlist)
                    out.write(output + '\n')
                elif name in pseudolist:
                    outlist = [seq, pval, name, 'pseudo']#, loc]
                    output = '\t'.join(outlist)
                    out.write(output + '\n')
                else:
                    outlist = [seq, pval, name, 'NA']#, loc]
                    output = '\t'.join(outlist)
                    o2.write(output + '\n')


