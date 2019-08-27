import os

'''
The alignments promoter names are compared to the ctss 
promoter name derived from the potential chromosomal location from 
the bamfile.
'''

####################################################################
def find_files(path, string):
    """Search a directory for files with string."""
    import os
    files = []
    dirs = os.listdir(path)
    for fle in dirs:
        if string in fle:
            files.append(fle)
    return files

####################################################################
def cross_reference_gene(infile):
    outfile = (infile.split('.'))[0] + '.' + (infile.split('.'))[1] + '.10.CR.tsv'
    o = open(outfile, 'w')
    with open(infile, 'r') as inF:
        for line in inF:
            if '@HD' in line:
                pass
            elif '@SQ' in line:
                pass
            elif '@PG' in line:
                pass
            else:
                linea = line.split('\t')
                bam_info = linea[0].split('|')
                bam_name = bam_info[6]
                bowtie_name = linea[2]

                if bam_name == bowtie_name:
                    o.write(line)

current = os.getcwd()
files = find_files(current, '.5.sam')
for infile in files:
    cross_reference_gene(infile)
