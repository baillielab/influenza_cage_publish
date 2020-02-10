
#leader_search_pilot_2016

from collections import Counter
import os
import time
import re
import random
import logging
from functools import reduce
#from slacker import Slacker
import multiprocessing as mp

path = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/10_other_datasets/'

file_list = []
dirs = os.listdir(path)
for fle in dirs:
    if 'genomicc_pilot_leader_list' in fle:
        file_list.append(fle)

print (file_list)

for infile in file_list:
	output = infile.split('.')
	leader_output = output[0] + '_tenmer.txt'
	background_list = []
	leader_list = []

	total = []
	with open(infile, 'r') as inF:
		sequence = []
		for line in inF:
			line = line.split('\t')
			if len(line[0]) >= 10:
				leader_tenmer = (line[0])[:10]
			else:
				pass
			sequence.append(leader_tenmer)
			leader_value = (line[1]).strip()
			leader_list.append([leader_tenmer, leader_value])
			total.append(int(leader_value))
	sequence_list = list(set(sequence))
	print (len(sequence_list))

	total_sequence = []
	for sequence in sequence_list:
		value = []
		for entry in leader_list:
			if sequence in entry:
				value.append(float(entry[1]))
		total_individual = sum(value)
		total_sequence.append([sequence, str(total_individual)])

	total_snatched = sum(total)

	ol = open(leader_output, 'w')
	for leader_tenmer in total_sequence:
		proportion = float(leader_tenmer[1]) / float(total_snatched)
		proportion = str(proportion * 100)
		ol.write(leader_tenmer[0] + '\t' + leader_tenmer[1] + '\t' + proportion + '\n')
