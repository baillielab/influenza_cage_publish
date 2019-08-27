import re
import os

def FluPromSearch(inputfile, prom, newpath):
#inputfile is in .sam format - the raw CAGE data
#prom is the promoter sequences
#newpath is the path of the output directory

	
	#get sample names
	in_file = os.path.abspath(inputfile)
	donor =  (re.split('%20|%3a|%29', inputfile))[8]
	sample =  (re.split('%20|%3a|%29', inputfile))[10]

	#new filename 
	o = (str(newpath) + '/%s_%s_%s.sam')%(prom,sample,donor)

	#check the file doesn't already exist
	if os.path.isfile(o) == False:

		#if it doesn't, open it
		o = open(o, 'w')

		#open the file and search for the promoter - not the special case for the SpliceProm is included below
		with open(in_file, 'r') as inF:
			for line in inF: 
				if prom in line:
					if prom == 'GCAAAAGCAG':
						if 'GCAAAAGCAGG' not in line:
							o.write(line)
					else:
						o.write(line)
	#if it does, print it and exit - that needs to be fixed
	else:
		print (('The file %s_%s_%s.sam already exists in the directory %s - the script will exit - please remedy')%(prom,sample,donor, str(newpath)))
		sys.exit()
