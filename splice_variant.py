##########################################################################################################################################################################################
## Created by: Sara Clohisey

'''
Orphans are influenza segment sequences (promoter adjacent on the 3' side)
that do not align at the correct position (12). 
There is an over-abundance of these sequences in the position range 9-15 
which are likely to be point mutations (polymerase) or sequencing error (HeliScope).
Those that align elsewhere may be true mRNA, with a leader and promoter, representing splice variants.
'''
## For Output:
import time
##dd/mm/yyyy format
date = (time.strftime("%Y%m%d"))


##########################################################################################################################################################################################
import re
import os
import regex
import sys
from collections import Counter
import numpy as np
from itertools import izip as zip, count



## Graphing

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import matplotlib.ticker as ticker

colors = ['#FF0000', '#FF8000', '#FFFF00', '#80FF00', '#00FFBF', '#00BFFF', '#0040FF', '#BF00FF', '#808080', 'k']
colors10 = ['k', '#FF0000', '#FF8000', '#FFFF00', '#80FF00', '#00FFBF', '#00BFFF', '#0040FF', '#BF00FF', '#808080', '#151B54'] 
greys = ['#d3d3d3', '#a8a8a8', '#7e7e7e', '#545454' ]
#plt.style.use('ggplot')

#Functions	##########################################################################################################################################################################################

def DictCreateFlu(libfile):
	l = open(libfile, 'r')
	linesl = l.readlines()
	l.close()
	lengl = len(linesl)
	dictionary = {}
	i = 0
	while i < lengl:
		if '>' in linesl[i]:
			seqid = str(linesl[i])
			if '>' not in linesl[i + 1]:
				seq = linesl[i + 1]
				seq = seq.replace('\n', '')
				seq = seq.replace('\r', '')
				dictionary[seqid] = seq
		else:
			pass
		i = i + 1
	return dictionary


def RegexSearchDictFluBasic(seq, dictionary):
	positions = []
	for i in range(len(dictionary.values())):
		m = re.search(seq, dictionary.values()[i]) 
		if m:
			positions.append([m.start(), dictionary.keys()[i]])

	return positions

#GetYourOrphans ###########################################################################################################################################################################################
segmented = '20170718_GCAAAGCAGG_segmented.txt'

orphans = []
with open(segmented, 'r',) as inF:
	for line in inF:
		if (line.split('\t'))[7] == 'Orphan':
			orphans.append((line.split('\t'))[5])

#Get the Genome as a dictionary###########################################################################################################################################################################################
libfile = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/resource/udorn_mrna_parsed.fa'

genome = DictCreateFlu(libfile)

counted = Counter(orphans)
counted = counted.most_common()

#Get Orphans which occur over x times ###########################################################################################################################################################################################
thresh = 10
orphans = []
orphans_count = []
i = 0
while i < len(counted):
	if int(counted[i][1]) > thresh:
		orphans.append(counted[i][0])
		orphans_count.append([counted[i][0], counted[i][1]])
	i = i + 1

#Search Genome for Orphans and output ###########################################################################################################################################################################################

outfile = (('orphan_possible_positions_%s.txt')%(thresh))

o = open(outfile, 'w')
o.write('Sequence\tCount\tSegment\tPosition\tMultimapGenome\n')
i = 0
multimaps = []
while i <len(orphans_count):
	seq = (orphans_count)[i][0]
	cnt = (orphans_count)[i][1]
	positions = RegexSearchDictFluBasic(seq, genome)
#Find Multimap ###########################################################################################################################################################################################
	for p in positions:
		o.write(('%s\t%s\t%s\t%s\t%s\n')%(seq, cnt, (p[1].replace('>', '')).replace('\n', ''), p[0], len(positions)))
	multimaps.append([seq, cnt, len(positions)])
	i = i + 1
o.close()

#Get Segments For Orphans ###########################################################################################################################################################################################
segments = []
with open(outfile, 'r') as inF:
	for line in inF:
		segments.append((line.split('\t'))[2])

segments = sorted(list(set(segments)))
segments =  segments[1:]

len_segments = len(segments)

#Graph For Orphans ###########################################################################################################################################################################################
tick_spacing_x = 500
tick_spacing_y = 25000
fig = plt.figure()
i = 0
while i < len_segments:
	segment = segments[i]
	length_key = len(genome['>' + segment + '\n'])
	
	x = []
	y = []
	s = []
	with open(outfile, 'r') as inF:
		for line in inF:
			if segment in line:
				x.append((line.split('\t'))[3])
				y.append((line.split('\t'))[1])
				s.append((10 - (int((line.split('\t'))[4])))**2)

	ax1 = fig.add_subplot(4,2,(i + 1)) #rows,cols,place
	ax1.scatter(x,y, s = s, color = colors[i])#, marker = '.')
	#ax1.scatter(x,y, color = colors[i], marker = '.')
	#ax1.tick_params(
    #axis='x',          # changes apply to the x-axis
    #which='both',      # both major and minor ticks are affected
    #bottom='off',      # ticks along the bottom edge are off
    #top='off',         # ticks along the top edge are off
    #labelbottom='off') # labels along the bottom edge are off
	
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_x))
	ax1.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
	ax1.set_title(segment, fontsize = 8)
	ax1.set_xlim(xmin = 14, xmax = length_key)
	ax1.set_ylim(ymin =5000, ymax =50000)

	i = i + 1
plt.tight_layout()
plt.show()
plt.savefig(('orphan_possible_positions_%s.pdf')%(thresh))

#Multimap to other orphans##########################################################################################################################################################################################

outfile2 = (('orphan_multimap_possible_positions_%s.txt')%(thresh))

o2 = open(outfile2, 'w')
o2.write('Sequence\tCount\tSegment\tPosition\tMultimapGenome\tMultimapOrphans\tAdjustedCount\n')

print count

for i in range(0, len(orphans)):
	value = list(map(lambda x: str(x).startswith(orphans[i]), (orphans)))
	indexes = ([a for a, j in zip(count(), value) if j == True])
	mm = list(map(lambda x : str((multimaps[x])[0]) + '\t' + str((multimaps[x])[1]) + '\t' + str((multimaps[x])[2]) , indexes))
	total = str(len(mm))
	with open(outfile, 'r') as inF:
		for line in inF:
			if orphans[i] == (line.split('\t'))[0]:
				adjust = float((line.split('\t'))[1]) / float(total)
				
				line = line.replace('\n', '')
				o2.write(line + '\t' + str(total) + '\t' +  str(adjust) + '\n')
		orphans_count[i].append(adjust)
o2.close()

print (orphans_count[0])[2]

#Add up adjusted multimap counts###########################################################################################################################################################################################

outfile3 = (('orphan_multimap_adjusttotalcount_possible_positions_%s.txt')%(thresh))

o3 = open(outfile3, 'w')
o3.write('Sequence\tCount\tSegment\tPosition\tMultimapGenome\tMultimapOrphans\tAdjustedCount\tTotalCountAdjusted\n')

for i in range(0, len(orphans)):
	value = list(map(lambda x: str(x).startswith(orphans[i]), (orphans)))
	indexes = ([a for a, j in zip(count(), value) if j == True])
	mm = list(map(lambda x: (orphans_count[x])[2], indexes))
	print mm
	mm = map(int, mm)
	total = sum(mm)
	with open(outfile2, 'r') as inF:
		for line in inF:
			if orphans[i] == (line.split('\t'))[0]:
				line = line.replace('\n', '')
				o3.write(line + '\t' + str(total) + '\n')
o3.close()

#Graph For Orphans ###########################################################################################################################################################################################
tick_spacing_x = 500
tick_spacing_y = 25000
fig = plt.figure()
i = 0
while i < len_segments:
	segment = segments[i]
	length_key = len(genome['>' + segment + '\n'])
	
	x = []
	y = []
	s = []
	with open(outfile3, 'r') as inF:
		for line in inF:
			if segment in line:
				x.append((line.split('\t'))[3])
				y.append((line.split('\t'))[7])
				s.append((line.split('\t'))[1])

	ax1 = fig.add_subplot(4,2,(i + 1)) #rows,cols,place
	ax1.scatter(x,y, s = s, color = colors[i])#, marker = '.')
	#ax1.scatter(x,y, color = colors[i], marker = '.')
	#ax1.tick_params(
    #axis='x',          # changes apply to the x-axis
    #which='both',      # both major and minor ticks are affected
    #bottom='off',      # ticks along the bottom edge are off
    #top='off',         # ticks along the top edge are off
    #labelbottom='off') # labels along the bottom edge are off
	
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_x))
	ax1.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing_y))
	ax1.set_title(segment, fontsize = 8)
	ax1.set_xlim(xmin = 14, xmax = length_key)
	#ax1.set_ylim(ymin =5000, ymax =50000)

	i = i + 1
plt.tight_layout()
plt.show()
plt.savefig(('orphan_possible_positions_multimapped_%s.pdf')%(thresh))
