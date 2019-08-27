#create bed file
import os

'''
Bedfile is created that states each 10mer
and each potential chromosomal
location for that tenmer.
'''

def create_bed(infile):
    current = os.path.dirname(os.path.realpath(__file__))
    infile_root = (infile.split('/'))[-1]
    outfile = '/' + (infile_root.split('.'))[0] + '.bed'
    o = open(current + outfile, 'w')

    abck_up_outfile = '/' + (infile_root.split('.'))[0] + '.filtered_out.bed'
    o2 = open(current + abck_up_outfile, 'w')


    with open(infile, 'r') as inF:
        next(inF)
        for line in inF:
            if 'sequence' in line:
                pass
            else:
                linea = (line.strip()).split('\t')
                if len(linea) >= 4:
                    tenmer = linea[0]
                    for location in linea[4:]:
                        for loc in location.split(", "):
                        #location = location.split(", ")
                            chrom = ((loc.split(' '))[0]).replace("['", '')
                            chrom = chrom.replace("'", '')
                            start = (((loc.split(' '))[1]).split(','))[0]
                            end = str(int(start) + 1)
                            strand = (((loc.split(','))[1]).split(':'))[0]
                            score = (((loc.split(','))[1]).split(':'))[1]
                            score = ((score.replace("'", '')).replace('[', '')).replace('\\n]', '')
                            outlist = [chrom, start, end, tenmer, score, strand]
                            output = '\t'.join(outlist)
                            o.write(output)
                            o.write('\n')
                else:
                    o2.write(line)


location_file = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/4_leaders/leader_tables/1a_chrom_locations/location_table_full.tsv'


create_bed(location_file)