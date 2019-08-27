#Header

##########################################################################################################################################################################################

#Purpose
'''
Conventional aligners can be specific enough for this process, 
however, I wanted something that was a bit more malleable. 
Also it doesn't take long so it was worth writing for the craic. 

The functions below allow for some variation in how the segment sequences, 
that is: sequences which follow the promoter sequence (GCAAAAGCAGG) 
in influenza mRNA (identified by promoter sequence), can be aligned 
to the influenza genome.

Exact searches are possible, only assigning a segment sequence to a segment 
if it aligns exactly, at the correct position, in this case 12 as the A 
which conventionally preceeds the promoter sequence has been excluded.
Additionally, the more verbose function allows for the user to have a 
more comprehensive output. It assigns exact matches, matches with a single mismatch.
In both cases if a sequence is unassigned it is classed as an 'orphan' 
and can be further processed in the variant search folder. 


Graphs can be made of the resulting segments so that segment dynamics over the timecourse can be observed.
'''

##########################################################################################################################################################################################
## Created by: Sara Clohisey
## Creation Date: 
## Last modified: 

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


## Graphing

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches

colors = ['#FF0000', '#FF8000', '#FFFF00', '#80FF00', '#00FFBF', '#00BFFF', '#0040FF', '#BF00FF', '#808080', 'k']
colors10 = ['k', '#FF0000', '#FF8000', '#FFFF00', '#80FF00', '#00FFBF', '#00BFFF', '#0040FF', '#BF00FF', '#808080', '#151B54'] 
greys = ['#d3d3d3', '#a8a8a8', '#7e7e7e', '#545454' ]
plt.style.use('ggplot')

###########################################################################################################################################################################################

#############################################CODE GOES HERE ##############################################################################################################################
#creates a dictionary with the libfile. This can then be done with the reverse file
#spliceprom and prom work with libfile while reverse promoter wrks with revlibfile
#runs inside RegexSearchDictFlu
def DictCreateFlu(libfile):


	#openfile and read in
	l = open(libfile, 'r')
	linesl = l.readlines()
	l.close()
	lengl = len(linesl)

	#make a dictionary ...which I apparently have named dictionary....
	dictionary = {}
	i = 0
		######### to print out if desired so you can check it later unindent and uncomment these.
		#o = 'udorngenomedict.txt'
		#o = open(o, 'w')
	while i < lengl:
		if '>' in linesl[i]:
			seqid = str(linesl[i])
			print seqid
			#o.write(seqid)
			if '>' not in linesl[i + 1]:
				seq = linesl[i + 1]
				seq = seq.replace('\n', '')
				seq = seq.replace('\r', '')
				print len(str([seq]))
				dictionary[seqid] = seq
				#o.write(seq)
		else:
			pass
		i = i + 1

		#print dictionary
		#o.close()
	return dictionary
#quick and messy to get the name from the current file
#not always needed
#run inside RegexSearchDictFlu

def FixNameFlu(seqfile):

	in_file = os.path.abspath(seqfile)
	in_file = in_file.split('.')
	in_file = in_file[0]
	in_file = in_file.split('/')
	print in_file
	name = (in_file[8].split('.'))[0]

	#print name
	return name

#seqfile is the promoter file
#name is from FixNameFlu
#dictionary is from DictCreateFlu
#EXAMPLE USE: RegexSearchDictFluBasic(file, FixNameFlu(file), DictCreateFlu(libfile))

#Basic is just a straightforward 'do these sequences align just after the promoter' question
#failing that the sequence is assigned as an orphan.
def RegexSearchDictFluBasic(seqfile, dictionary):

	import re
	import os

	with open(seqfile, 'r') as inF:
		#LINE: [0]Promoter 	[1]	TAG# 	[2]SampleName 	[3]SequenceTAG	[4]	SequenceLen 	[5]LeaderSeq		[6]LeaderLen 	[7]SegmentSeq 	[8]SegmentLen 		[9]QualityScore 	[10]Donor 	[11]Hour 	[12]Treatment

		o = (('%s_%s_segmented.txt')%((str(date), ((str(seqfile).split('.'))[-1]))))
		o = open(o, 'w')
		o.write('Promoter\tTAG#\tSampleName\tLeader\tLeaderLen\tSegmentSeq\tSegmentLen\tSegment#\tSegmentPostn\n')

		for line in inF:
			#skips header if CAGE data
			if 'Promoter' not in line: 
				#define variable for each iteration
				line = line.split('\t')
				segseq = str(line[7])
				#sequences with less than three nts are unassignable - in fact should be removed in new code.
				if len(segseq) < 3:
					o.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(line[0],line[1],line[2],line[5],line[6],line[7],line[8],'Short','X'))
					#pass
				else:
					temp = {}

					#search the dictionary for the segsequence
					for i in range(len(dictionary.values())):
						m = re.search(segseq, dictionary.values()[i]) 
						if m:
							#print m.start()
							temp[m.start()] = dictionary.keys()[i]
					#print temp

					#the following takes into account that position 12 is the most preferable position. I tried to do this from a list but ran into some problems that were dealt with easier in the messy code...

					if temp.has_key(12):
						value = (str(temp[12]).replace('\n', '')).replace('>', '')
						o.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(line[0],line[1],line[2],line[5],line[6],line[7],line[8],str(value),12))

					#if the segseq does not fit with any of these it becomes an orphan. Orphans must be later aligned (smaligned) to determine which are most likely
				
					else:
						o.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(line[0],line[1],line[2],line[5],line[6],line[7],line[8],'Orphan','X'))

		o.close()

#seqfile is the promoter file
#name is from FixNameFlu
#dictionary is from DictCreateFlu
#EXAMPLE USE: RegexSearchDictFluVerbose(file, FixNameFlu(file), DictCreateFlu(libfile))
#file = (('collate_%s.txt')%(promoter))

#Verbose is a more complicated question. If the exact match is not found, the function then checks if a single mismatch can lead to an alignment
#if that is not the case the function checks if the segment sequence contains part of the promoter sequence.
#failing that the sequence is assigned as an orphan.
def RegexSearchDictFluVerbose(seqfile, dictionary):

	import re
	import os
	import regex


	prom = [0,1,2,3,4,5,6,7,8,9,10,11]
	with open(seqfile, 'r') as inF:
		#LINE: [0]Promoter 	[1]	TAG# 	[2]SampleName 	[3]SequenceTAG	[4]	SequenceLen 	[5]LeaderSeq		[6]LeaderLen 	[7]SegmentSeq 	[8]SegmentLen 		[9]QualityScore 	[10]Donor 	[11]Hour 	[12]Treatment

		o = o = (('%s_%s_segmented_verbose.txt')%((str(date), ((str(seqfile).split('.'))[-1]))))
		o = open(o, 'w')
		o.write('Promoter\tTAG#\tSampleName\tLeader\tLeaderLen\tSegmentSeq\tSegmentLen\tSegment#\tSegmentPostn\n')

		for line in inF:
			#SKIPS HEADER IF CAGE DATA
			if 'Promoter' not in line: 
				#DEFINE VARIABLE FOR EACH ITERATION
				line = line.split('\t')
				segseq = str(line[7])
				#SEQUENCES WITH LESS THAN THREE NTS ARE UNASSIGNABLE TO A 
				if len(segseq) < 3:
					o.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(line[0],line[1],line[2],line[5],line[6],line[7],line[8],'Short','X'))
				else:
					#CREATE A TEMPORARY DICTIONARY
					temp = {}
					#SEARCH THE DICTIONARY FOR THE SEGSEQUENCE
					for i in range(len(dictionary.values())):
						m = re.search(segseq, dictionary.values()[i]) 
						if m:
							temp[m.start()] = dictionary.keys()[i]

					#THE FOLLOWING TAKES INTO ACCOUNT THAT POSITION 12 IS THE MOST PREFERABLE POSITION.
					#ASSIGN PERFECTLY ALIGNED SEQUENCES TO THE FIRST THREE BASES OF THE SEGMENT

					if temp.has_key(12): 
						value = (str(temp[12]).replace('\n', '')).replace('>', '')
						o.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(line[0],line[1],line[2],line[5],line[6],line[7],line[8],str(value),12))
					#IF CAN'T PERFECTLY ALIGN TO SEGMENT START
					#MAKE USE OF MISMATCHING AND SEE IF THE SEQUENCE MATCHES WITH A SINGLE MISMATHC. THIS CODE DOES NOT DIFFERENTIATE BETWEEN A SUB, INSERT OR DEL. IT'S SIMPLY A 1 BASE MISMATCH.

					else:
						for i in range(len(dictionary.values())):
							pattern = (('(%s){e<=1}')%(segseq))
							m = regex.search(pattern, dictionary.values()[i]) 
							if m:
								temp[m.start()] = dictionary.keys()[i]
						#ASSIGN 1MM ALIGNED SEQUENCES TO THE FIRST THREE BASES OF THE SEGMENT
						
						if temp.has_key(12):
							value = (str(temp[12]).replace('\n', '')).replace('>', '')
							o.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(line[0],line[1],line[2],line[5],line[6],line[7],line[8],str(value),12))

						#IF THE SEQUENCE DOES NOT ALIGN IN THESE CASES, CHECK IF IT IS A CASE OF 'PRIME AND REALIGN' IE THAT THE SEGMENT SEQUENCE CONTAINS PART OF THE PROMOTER.
					#else:
					#	for p in prom:
					#		if temp.has_key(int(p)):
					#			value = (str(temp[int(p)].replace('\n', '')).replace('>', ''))
					#			o.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(line[0],line[1],line[2],line[5],line[6],line[7],line[8],'PAR',str(p)))
								#break

						#FINALLY IF  THE SEGSEQ DOES NOT FIT WITH ANY OF THESE IT BECOMES AN ORPHAN. ORPHANS MUST BE LATER ALIGNED TO DETERMINE WHICH ARE MOST LIKELY TO BE TRUE SEQUENCES
				
						else:
							o.write(('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n')%(line[0],line[1],line[2],line[5],line[6],line[7],line[8],'Orphan','X'))

###########################################################################################################################################################################################
#Graphs

#uses output from RegexSearchDictFluExactSeg
#LINE: [0]Promoter 	[1]TAG 	[2SampleName 	[3]Leader 	[4]LeaderLen 	[5]SegmentSeq 	[6]SegmentLen 	[7]Segment# 	[8]SegmentPostn'
#creates bar charts
#can be used to graph all at once or one at a time
def GraphSegmentsHour(seqfile):
	from collections import Counter
	import numpy as np
	from itertools import chain


	#barchart of segments by hour

	hourlst = ['0','2','7','24']

	seglist = ['Orphan', 'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5', 'Segment 6', 'Segment 7', 'Segment 8', 'Short']

	segarraylist = []
	percentarray = []
	for h in hourlst:
		print h
		segments = []

		with open(seqfile, 'r') as inF:
			for line in inF:	
				if 'Promoter' not in line:
					line = line.split('\t')
					
					
					if h == ((line[2].split('_'))[2]).replace('h', ''):
						segment = line[7]
						segments.append(segment)

		countedseg = Counter(segments)
		countedsegordered = sorted(countedseg.items())



		segs = []
		i = 0
		while i < len(seglist):

			segs.append(int((((str(countedsegordered[i])).split(','))[1]).replace(')', '')))
			i = i + 1

		my_order = [1,2,3,4,5,6,7,8,9,0]
		segs = [segs[i] for i in my_order]


		segarraylist.append(segs)
		segaxis= ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan')
		#plot 

		percentseg = sum(segs)



		percentseglist = []

		for i in segs:
			i = (float(i) / float(percentseg)) * 100

			percentseglist.append(i)

		ass = sum(percentseglist[0:8])

		unass = sum(percentseglist[8:])

		if float(ass) + float(unass) != 100:
			break

		percentarray.append(percentseglist)

		plt.clf()
		#input
		y = segs
		
		x_pos = np.arange(len(segaxis))
		#labels
		name = (('Segment_%shours_Raw')%(h))	
		plt.ylabel('Reads')
		plt.xlabel('Segment')
		plt.title(('Segments at %s hours')%(h))
		#layout
		plt.grid(True)
		width = 0.75	
		plt.bar(x_pos, y, width = width, align='center')
		plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
		plt.xlim(xmin=-1)
		plt.tight_layout()
		#show
		plt.show()
		#export as pdf
		plt.savefig(name)

		plt.clf()
		#input
		y = percentseglist
		
		x_pos = np.arange(len(segaxis))
		#labels
		name = (('Segment_%shours_Percent')%(h))	
		plt.ylabel('Percent (%)')
		plt.xlabel('Segment')
		plt.title(('Percent Segments at %s hours')%(h))
		#layout
		plt.grid(True)
		width = 0.75	
		plt.bar(x_pos, y, width = width, align='center')
		plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
		plt.xlim(xmin=-1)
		plt.tight_layout()
		#show
		plt.show()
		#export as pdf
		plt.savefig(name)






	#colors = ['r','g','b','y']

	print percentarray

	a = np.array(segarraylist)
	b = np.array(percentarray)

	print b
	#print a
	#print a.ndim
	#print a.mean(axis=0) 

	y = a.mean(axis=0)
	std = np.std(a, axis=0)
	#print std
	name = ('/segments/Average_Segment_hours')	
	plt.ylabel('Average_Reads')
	plt.xlabel('Segment')
	plt.title('Average Segments Overall')
	#layout
	plt.grid(True)
	width = 0.75	
	plt.bar(x_pos, y, width = width, align='center', yerr=std)
	plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
	plt.xlim(xmin=-1)
	plt.tight_layout()
	#show
	plt.show()
	#export as pdf
	plt.savefig(name)


	plt.clf()

	y = b.mean(axis=0)

	std = np.std(b, axis=0)
	#print std
	name = ('Average_Segment_Percent_hours')	
	plt.ylabel('Percent (%)')
	plt.xlabel('Segment')
	plt.title('Average Percent Segments Overall')
	#layout
	plt.grid(True)
	width = 0.75	
	plt.bar(x_pos, y, width = width, align='center', yerr=std)
	plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
	plt.xlim(xmin=-1)
	plt.tight_layout()
	#show
	plt.show()
	#export as pdf
	plt.savefig(name)

#uses output from RegexSearchDictFluExactSeg
#LINE: [0]Promoter 	[1]TAG 	[2SampleName 	[3]Leader 	[4]LeaderLen 	[5]SegmentSeq 	[6]SegmentLen 	[7]Segment# 	[8]SegmentPostn'
#creates bar charts
#barchart of segments by hour
#can be used to graph all at once or one at a time
def GraphSegmentsDonor(seqfile):

	from collections import Counter
	import numpy as np
	from itertools import chain
	

	donorlst = ['1', '2', '3', '4']

	seglist = ['Orphan', 'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5', 'Segment 6', 'Segment 7', 'Segment 8', 'Short']

	segarraylist = []
	percentarray = []
	for d in donorlst:
		print d
		segments = []

		with open(seqfile, 'r') as inF:
			for line in inF:	
				if 'Promoter' not in line:
					line = line.split('\t')
					
					if d == (line[2])[-1]:
						segment = line[7]
						segments.append(segment)

		countedseg = Counter(segments)
		countedsegordered = sorted(countedseg.items())



		segs = []
		i = 0
		while i < len(seglist):
			segs.append(int((((str(countedsegordered[i])).split(','))[1]).replace(')', '')))
			i = i + 1

		my_order = [1,2,3,4,5,6,7,8,9,0]
		segs = [segs[i] for i in my_order]


		segarraylist.append(segs)
		segaxis= ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan')
		#plot 

		percentseg = sum(segs)



		percentseglist = []

		for i in segs:
			i = (float(i) / float(percentseg)) * 100

			percentseglist.append(i)

		ass = sum(percentseglist[0:8])

		unass = sum(percentseglist[8:])

		if float(ass) + float(unass) != 100:
			break

		percentarray.append(percentseglist)

		plt.clf()
		#input
		y = segs
		
		x_pos = np.arange(len(segaxis))
		#labels
		name = (('Segment_donor_%s_Raw')%(d))	
		plt.ylabel('Reads')
		plt.xlabel('Segment')
		plt.title(('Segment Distribution Donor %s')%(d))
		#layout
		plt.grid(True)
		width = 0.75	
		plt.bar(x_pos, y, width = width, align='center')
		plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
		plt.xlim(xmin=-1)
		plt.tight_layout()
		#show
		plt.show()
		#export as pdf
		plt.savefig(name)

		plt.clf()
		#input
		y = percentseglist
		
		x_pos = np.arange(len(segaxis))
		#labels
		name = (('Segment_donor_%s_Percent')%(d))	
		plt.ylabel('Percent (%)')
		plt.xlabel('Segment')
		plt.title(('Percent Segment Distribution Donor %s')%(d))
		#layout
		plt.grid(True)
		width = 0.75	
		plt.bar(x_pos, y, width = width, align='center')
		plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
		plt.xlim(xmin=-1)
		plt.tight_layout()
		#show
		plt.show()
		#export as pdf
		plt.savefig(name)

	a = np.array(segarraylist)
	b = np.array(percentarray)



	y = a.mean(axis=0)
	std = np.std(a, axis=0)
	#print std
	name = ('Average_Segment_Donor')	
	plt.ylabel('Average_Reads')
	plt.xlabel('Segment')
	plt.title('Average Segments Overall')
	#layout
	plt.grid(True)
	width = 0.75	
	plt.bar(x_pos, y, width = width, align='center', yerr=std)
	plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
	plt.xlim(xmin=-1)
	plt.tight_layout()
	#show
	plt.show()
	#export as pdf
	plt.savefig(name)


	plt.clf()

	y = b.mean(axis=0)





	std = np.std(b, axis=0)

	name = ('Average_Segment_Percent_Donor')	
	plt.ylabel('Percent (%)')
	plt.xlabel('Segment')
	plt.title('Average Percent Segments Overall')
	#layout
	plt.grid(True)
	width = 0.75	
	plt.bar(x_pos, y, width = width, align='center', yerr=std)
	plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
	plt.xlim(xmin=-1)
	plt.tight_layout()
	#show
	plt.show()
	#export as pdf
	plt.savefig(name)

#uses output from RegexSearchDictFluExactSeg
#LINE: [0]Promoter 	[1]TAG 	[2SampleName 	[3]Leader 	[4]LeaderLen 	[5]SegmentSeq 	[6]SegmentLen 	[7]Segment# 	[8]SegmentPostn'
#creates bar charts
#can be used to graph all at once or one at a time
def GraphSegmentsDonorHour(seqfile):


	from collections import Counter
	import numpy as np
	from itertools import chain

	donorlst = ['1', '2', '3', '4']

	hourlst = ['0', '2', '7', '24']

	seglist = ['Orphan', 'Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5', 'Segment 6', 'Segment 7', 'Segment 8', 'Short']

	segarraylist = []
	percentarray = []
	for d in donorlst:
		for h in hourlst:
			print d
			print h
			segments = []

			with open(seqfile, 'r') as inF:
				for line in inF:	
					if 'Promoter' not in line:
						line = line.split('\t')
						
						if d == (line[2])[-1] and h == ((line[2].split('_'))[2]).replace('h', ''):
							segment = line[7]
							segments.append(segment)

			countedseg = Counter(segments)
			countedsegordered = sorted(countedseg.items())



			segs = []
			i = 0
			while i < len(seglist):
				segs.append(int((((str(countedsegordered[i])).split(','))[1]).replace(')', '')))
				i = i + 1

			my_order = [1,2,3,4,5,6,7,8,9,0]
			segs = [segs[i] for i in my_order]


			segarraylist.append(segs)
			segaxis= ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan')
			#plot 

			percentseg = sum(segs)



			percentseglist = []

			for i in segs:
				i = (float(i) / float(percentseg)) * 100

				percentseglist.append(i)

			ass = sum(percentseglist[0:8])

			unass = sum(percentseglist[8:])

			if float(ass) + float(unass) != 100:
				break

			percentarray.append(percentseglist)

			plt.clf()
			#input
			y = segs
			
			x_pos = np.arange(len(segaxis))
			#labels
			name = (('Segment_donor_%s_hour_%s_Raw')%(d,h))	
			plt.ylabel('Reads')
			plt.xlabel('Segment')
			plt.title(('Segment Distribution Donor %s Hour %s')%(d,h))
			#layout
			plt.grid(True)
			width = 0.75	
			plt.bar(x_pos, y, width = width, align='center')
			plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
			plt.xlim(xmin=-1)
			plt.tight_layout()
			#show
			plt.show()
			#export as pdf
			plt.savefig(name)

			plt.clf()
			#input
			y = percentseglist
			
			x_pos = np.arange(len(segaxis))
			#labels
			name = (('Segment_donor_%s_hour_%s_Percent')%(d,h))	
			plt.ylabel('Percent (%)')
			plt.xlabel('Segment')
			plt.title(('Percent Segment Distribution Donor %s Hour %s')%(d,h))
			#layout
			plt.grid(True)
			width = 0.75	
			plt.bar(x_pos, y, width = width, align='center')
			plt.xticks(x_pos, ('1', '2', '3', '4', '5','6','7','8', 'Short', 'Orphan'), rotation='vertical')
			plt.xlim(xmin=-1)
			plt.tight_layout()
			#show
			plt.show()
			#export as pdf
			plt.savefig(name)
###########################################################################################################################################################################################

#run code

#Segment the file
promoterfile = 'GCAAAAGCAGG_CAGE_Sequences.txt'
libfile = 'udorn_mrna_parsed.fa'

###########################################################################################################################################################################################


#Footer
###########################################################################################################################################################################################
