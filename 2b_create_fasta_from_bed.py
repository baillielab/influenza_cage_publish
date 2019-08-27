#create_fasta
import os
'''Converts overlap output to fasta for downstream analysis.'''


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
def create_fasta(infile):
    with open(infile, 'r') as inF:
        outfile = (infile.split('.'))[0] + '.' + (infile.split('.'))[1] + '.fasta'
        o = open(outfile, 'w')
        for line in inF:
            linea = line.split('\t')
            print linea
            outline1 = '|'.join(['>', linea[0], linea[1], linea[2], linea[4], linea[5], linea[9], linea[10], '\n'])
            o.write(outline1)
            outline2 = linea[3] + '\n'
            o.write(outline2)

####################################################################
current = os.getcwd()
files = find_files(current, '.5.overlap')

for infile in files:
    create_fasta(infile)