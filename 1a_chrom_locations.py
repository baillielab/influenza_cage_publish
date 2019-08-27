#retrieve chrom locations
import os

def find_files(path, string):
    import os
    files = []
    dirs = os.listdir(path)
    for fle in dirs:
        if string in fle:
            files.append(fle)
    return files
capsnatch = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
leaders = '4_leaders/publish/'
tenmer_locations_file = 'rename_10mers_no_database/0_tenmer_locations/'
tenmer_locations_full = capsnatch + leaders + tenmer_locations_file
tenmer_table_file = 'leader_tables/1_tenmers/tenmer_table_full.tsv'
tenmer_table = capsnatch + '4_leaders/' + tenmer_table_file
#######################################################
tenmer_locations_files = sorted(find_files(tenmer_locations_full, '_tenmers_b.tsv'))
print tenmer_locations_files

outfile = 'location_table_full.tsv'
o = open(outfile, 'w')
o.write('tenmer\thour\tdonor\ttotal in sample\tlocations\n')

with open(tenmer_table, 'r') as inF:
    next(inF)
    for line in inF:
        linea = line.split('\t')
        tenmer = linea[0]
        donor = linea[2].replace('d','donor')
        hour = linea[1]
        name = (('Ud_%s_%s')%(hour, donor))
        for fle in tenmer_locations_files:
            if name in fle:
                fle = tenmer_locations_full + fle
                with open(fle, 'r') as inF2:
                    for line2 in inF2:
                        if tenmer in line2:
                            line2a = line2.split('\t')
                            o.write(tenmer + '\t' + hour + '\t' + donor + '\t' + line2a[1] + '\t' + str(line2a[2:]) + '\n')
