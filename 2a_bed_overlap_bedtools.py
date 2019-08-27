import os

'''
Bedfile is overlapped with hg19, with a 5bp window +/- and exact strand 
match only.
'''

infile = 'location_table_full.bed'
outfile = (infile.split('.'))[0] + '5.overlap'
cmd = (('bedtools window -a %s -b hg19.bed -w 5 -sm > %s')%(infile, outfile))
os.system(cmd)