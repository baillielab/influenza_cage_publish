
#leader_search_bartel

from collections import Counter
import os
import time
import re
import random
import logging
from functools import reduce
#from slacker import Slacker
import multiprocessing as mp

path = '/Users/sclohise/Desktop/'

def Sikora(file_list):

	for infile in file_list:
		output = infile.split('.')
		leader_output = output[0] + '_leader_tenmer.txt'
		background_output = output[0] + '_background_tenmer.txt'
		
		background_list = []
		flu_list = []


		with open(infile, 'r') as inF:
			for line in inF:
				if '>' in line:
					pass
				elif 'GCAAAAGCAGG' in line:
					flu_list.append(line[:10])
				else:					
					background_list.append(line[:10])

		print (background_list)

		o1 = open(leader_output, 'w')
		counted = Counter(flu_list)
		newv = Counter.most_common(counted)
		for x in newv:
			o1.write(x[0] + '\t' + str(x[1]) +  '\n')

		o = open(background_output, 'w')
		counted = Counter(background_list)
		newv = Counter.most_common(counted)
		for x in newv:
			o.write(x[0] + '\t' + str(x[1]) +  '\n')


def Gu(file_list):

	for infile in file_list:
		output = infile.split('.')
		leader_output = output[0] + '_leader_tenmer.txt'
		tenmer_output = output[0] + '_background_tenmer.txt'
		flu_output = output[0] + '_IVA.txt'
		background_list = []
		flu_list = []
		leader_list = []

		if 'IVA' in infile:
			with open(infile, 'r') as inF:
				for line in inF:
					if '#' in line:
						pass
					elif 'GCAAAAGCAGG' in line:
						flu_list.append(line.strip())
						line = line.split('\t')
						leader_tenmer = (line[0])[:10]
						leader_value = (line[1]).strip()
						leader_list.append([leader_tenmer, leader_value])

					else:
						tenmer = ((line.split('\t'))[0])[:10]
						value = ((line.split('\t'))[1]).strip()
						background_list.append([tenmer, value])

				o = open(flu_output, 'w+')
				for seq in flu_list:
					o.write(seq + '\n')

				ol = open(leader_output, 'w')
				for leader_tenmer in leader_list:
					ol.write(leader_tenmer[0] + '\t' + leader_tenmer[1] + '\n')

				o1 = open(tenmer_output, 'w')
				for tenmer in background_list:
					o1.write(tenmer[0] + '\t' + tenmer[1] + '\n')

		else:
			with open(infile, 'r') as inF:
				for line in inF:
					if '#' in line:
						pass
					else:
						tenmer = ((line.split('\t'))[0])[:10]
						value = ((line.split('\t'))[1]).strip()
						background_list.append([sequence, value])
			o1 = open(tenmer_output, 'w')
			for tenmer in background_list:
				o1.write(tenmer[0] + '\t' + tenmer[1] + '\n')


def Bartel(file_list):

	for infile in file_list:
		output = infile.split('.')
		leader_output = output[0] + 'leader.txt'
		tenmer_output = output[0] + 'tenmer.txt'
		leaders_list = []
		tenmers_list = []
		with open(infile, 'r') as inF:
			for line in inF:
				if 'ORIGINAL_SEQUENCE' in line:
					line_a = line.split('\t')
					place = line_a.index('ORIGINAL_SEQUENCE')
		with open(infile, 'r') as inF:
			for line in inF:
				leader = (line.split('\t'))[place]
				leaders_list.append(leader)

				if len(leader) >= 10:
					tenmers_list.append(leader[:10])

		o = open(leader_output, 'w+')
		o.write(str(len(leaders_list)) + '\n')
		for leader_seq in leaders_list:
			o.write(leader_seq + '\n')

		o1 = open(tenmer_output, 'w+')
		o1.write(str(len(tenmers_list)) + '\n')

		counted = Counter(tenmers_list)
		newv = Counter.most_common(counted)
		for x in newv:
			o1.write(x[0] + '\t' + str(x[1]) +  '\n')
	

def find_leaders_in_other_datasets(author):

	if author == 'gu':
		search = 'GSM1648'

	elif author == 'bartel'
		search = ''

	elif author == 'sikora':
		search = 'sra'

	else:
		print ('Wrong author - no dataset')

	file_list = []
	dirs = os.listdir(path)
	for fle in dirs:
	    if search in fle:
	        file_list.append(fle)

	print (file_list)


file_list = find_leaders_in_other_datasets(author)

