tenmer_table_full = {}
tenmer_table_stats = {}
#######################################################
capsnatch = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
leaders = '4_leaders/'
table = '/leader_tables/1_tenmers/tenmer_table_full_corrected.tsv'
#######################################################
tenmer_table_full_file = capsnatch + leaders + table
with open(tenmer_table_full_file, 'r') as inF:
    next(inF)
    for line in inF:
        linea = (line.strip()).split('\t')
        name = linea[0] + '_' + linea[1] + '_' + linea[2]
        #pvalue_OR
        tenmer_table_stats[name] = linea[3] + '_' + linea[4]
        #snatched_unsnatched
        tenmer_table_full[name] = linea[5] + '_' + linea[6]

tenmer_location = {}
#######################################################
location = 'leader_tables/1a_chrom_locations/location_table_full.tsv'
#######################################################
tenmer_location_full_file = capsnatch + leaders + location
with open(tenmer_location_full_file, 'r') as inF:
    next(inF)
    for line in inF:
        linea = (line.strip()).split('\t')
        if len(linea) > 2:
            name = linea[0] + '_' + linea[1] + '_' + linea[2].replace('onor', '')
            tenmer_location[name] = linea[3:]
        else:
            pass
#######################################################
print tenmer_table_full['TCGGAGCGGT_0h_d1']
print tenmer_table_stats['TCGGAGCGGT_0h_d1']
print tenmer_location['TCGGAGCGGT_0h_d1']

output_file = 'collated_tenmers_locations.tsv'
o = open(output_file, 'w')
o.write('tenmer\thour\tdonor\tFDR\tOR\tsnatched\tunsnacthed\ttotal\tlocations\n')
with open(tenmer_table_full_file, 'r') as inF:
    next(inF)
    for line in inF:
        linea = (line.strip()).split('\t')
        name = linea[0] + '_' + linea[1] + '_' + linea[2]
        stats = tenmer_table_stats[name]
        full = tenmer_table_full[name]
        try:
            loc = str(tenmer_location[name])
        except KeyError:
            loc = 'None'

        stats = stats.split('_')
        full = full.split('_')

        outlist = [linea[0], linea[1], linea[2], stats[0], stats[1], full[0], full[1], loc, '\n']
        output = '\t'.join(outlist)
        o.write(output)