import os
#choose a gene name

def find_files(path, string):
    import os
    files = []
    dirs = os.listdir(path)
    for fle in dirs:
        if string in fle:
            files.append(fle)
    return files

#######################################################
capsnatch = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
leaders = '4_leaders/publish/'
indir = 'fishers_unbiased_leader_search/uncollated_donors_repeat/adjusted_pvalues/'
thisseq_not_snatched = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/10mercount_tables/thisseq_notsnatched'
thisseq_snatched = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/10mercount_tables/thisseq_snatched'

#######################################################

list_of_tenmers = []
tenmer_count_tables = 'source_data/10mercount_tables/'
tenmer_count_tables_full = capsnatch + tenmer_count_tables
list_tenmers_file = tenmer_count_tables_full + 'thisseq_snatched'

with open(list_tenmers_file, 'r') as inF:
    next(inF)
    for line in inF:
        linea = line.split('\t')
        list_of_tenmers.append(linea[0])
#######################################################

#/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/all_tenmers
all_tenmers = capsnatch + '/source_data/all_tenmers/'
adjusted_pvalue = capsnatch + leaders + indir
tenmer_fishers = sorted(find_files(adjusted_pvalue, 'adjusted_pvalues.txt'))
SC_tenmer_counts = sorted(find_files(all_tenmers, '_tenmers_SC.tsv'))

'''
###########################################
outfile = 'tenmer_table_full.tsv'
o = open(outfile, 'w')
o.write('tenmer\thour\tdonor\tOR\tFDR\tsnatch\tnotsnatch\n')

for tenmer in list_of_tenmers:
    print tenmer
    i = 0
    while i < len(tenmer_fishers):
        hour = (tenmer_fishers[i].split('_'))[1]
        donor = ((tenmer_fishers[i].split('_'))[2]).replace('adjusted','')
        infile = all_tenmers + SC_tenmer_counts[i]
        with open(infile, 'r') as inF:
            for line in inF:
                if tenmer in line:
                    linea = (line.strip()).split('\t')
                    unsnatch = linea[1]
                    snatch = linea[2]
        infile = capsnatch + leaders + indir + tenmer_fishers[i]
        f = open(infile, 'r')
        lines = f.readlines()
        f.close()
        if tenmer in lines:
            for line in lines:
                if tenmer in line:
                    linea = line.split('\t')
                    FDR = linea[2]
                    OR = linea[3]
        else:
            FDR = 'NA'
            OR = 'NA'

        outlist = [tenmer, hour, donor, OR, FDR]
        outlist = outlist + [snatch, unsnatch, '\n']
        output = '\t'.join(outlist)
        o.write(output)
        i = i + 1
#######################################################
'''

tenmer_table_full = 'tenmer_table_full_test.tsv'
output = 'tenmer_table_full_corrected.tsv'
o = open(output, 'w')
o.write('tenmer\thour\tdonor\tOR\tFDR\tsnatch\tnotsnatch\n')

for fle in tenmer_fishers:
    name = fle.split('_')
    donor = (name[2]).replace('adjusted', '')
    hour = name[1]
    dictname = hour + donor
    index_number = tenmer_fishers.index(fle)
    index_number = (int(index_number) + 9)

    with open(thisseq_not_snatched, 'r') as inF:
        for line in inF:
            linea = (line.strip()).split('\t')
            snatched = linea[index_number]

    dictname = {}

    infile = capsnatch + leaders + indir + fle
    with open(infile, 'r') as inF:
        next(inF)
        next(inF)
        for line in inF:
            linea = line.split('\t')
            dictname[linea[0]] = linea[2] + '\t' + linea[3]
    with open(tenmer_table_full, 'r') as inF2:
        next(inF2)
        for line in inF2:
            if donor in line:
                if hour in line:
                    linea = line.split('\t')
                    tenmer = linea[0]
                    try:
                        output = dictname[linea[0]]
                    except KeyError:
                        output = 'NA\tNA'
                    print output

                    with open(thisseq_not_snatched, 'r') as inF:
                        for line in inF:
                            if tenmer in line:
                                linea = (line.strip()).split('\t')
                                notsnatched = linea[index_number]
                    with open(thisseq_snatched, 'r') as inF:
                        for line in inF:
                            if tenmer in line:
                                linea = (line.strip()).split('\t')
                                snatched = linea[index_number]
                    o.write(tenmer + '\t' + hour + '\t' + donor + '\t' + output + '\t' + snatched + '\t' + notsnatched + '\n')





'''
with open(infile, 'r') as inF:
    for line in inF:
        linea = line.split('\t')
'''

