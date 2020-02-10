###############################################################################
from itertools import izip as zip, count
from collections import Counter
import os
import time
import re
import random
import logging
from functools import reduce
import multiprocessing as mp
###############################################################################

###############################################################################
tenmer_path = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
file_name = 'tenmer_SC.txt'
tenmer_file = tenmer_path + file_name
tenmer_list = []

print tenmer_file

with open(tenmer_file, 'r') as inF:
	for line in inF:
		tenmer = line.strip()
		tenmer_list.append(tenmer)

#initialise file
output = 'other_dataset_summary_table.txt'
o = open(output, 'w')
o.write('tenmer' + '\n')
for tenmer in tenmer_list:
	o.write(tenmer + '\n')
o.close()

f = open('other_dataset_summary_table.txt', 'r')
lines = f.readlines()
f.close()

outlines = []
for line in lines:
	outlines.append([line])

###############################################################################
def filter_value( someList, value ):
    for x, y in someList:
        if x == value :
            yield x,y

###############################################################################
path = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/10_other_datasets'
file_list = []
dirs = os.listdir(path)
for fle in dirs:
    if '_tenmer.txt' in fle:
        file_list.append(fle)
file_list = sorted(file_list)

print file_list


for infile in file_list:
	print infile
	output_list = []
	with open(infile, 'r') as inF:
		for line in inF:
			line = (line.strip()).split('\t')
			if len(line) < 2:
				pass
			else:
				sequence = line[0]
				percent = line[2]
				output_list.append([sequence, percent])
	print len(output_list)
###############################################################################
	i = 0
	while i < len(lines):
		print i
		sequence = ((lines[i].strip()).split('\t'))[0]
		print sequence 
		result = list(filter_value(output_list, sequence))
		print result
		if len(result) >=1:
			out = [x[1] for x in result]
			(outlines[i]).append(out)			
		else:
			(outlines[i]).append('NA')
		i = i + 1


output = 'other_dataset_summary_table.txt'
o = open(output, 'w')
o.write('tenmer')
for name in file_list:
	o.write('\t' + name )
o.write('\n')

for out_line in outlines:
	output = '\t'.join(out_line)
	o.write(str(output))
	o.write('\n')
		



###############################################################################