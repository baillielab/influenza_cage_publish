from itertools import izip as zip, count
from collections import Counter
import os
import time
import re
import logging
from functools import reduce
import multiprocessing as mp


###############################################################################
# get locations
###############################################################################
output_path = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/'
path = '/mnt/ris-fas1a/linux_groups2/fantom5/Clustering/FLU_TIMECOURSE/bamfiles/'
###############################################################################
unfinished = ['fishers_output_Ud_2h_d1_filteredfdr.txt',
              'fishers_output_Ud_2h_d4_filteredfdr.txt',
              'fishers_output_Ud_7h_d2_filteredfdr.txt',
              'fishers_output_Ud_7h_d3_filteredfdr.txt']

finished = ['fishers_output_Ud_2h_d2_filteredfdr.txt',
            'fishers_output_Ud_2h_d3_filteredfdr.txt',
            'fishers_output_Ud_7h_d1_filteredfdr.txt',
            'fishers_output_Ud_7h_d4_filteredfdr.txt',
            'fishers_output_Ud_24h_d2_filteredfdr.txt',
            'fishers_output_Ud_24h_d3_filteredfdr.txt',
            'fishers_output_Ud_24h_d1_filteredfdr.txt',
            'fishers_output_Ud_24h_d4_filteredfdr.txt']

def get_files(path, string, string_2):
    files = []
    dirs = os.listdir(path)
    for fle in dirs:
        if string in fle:
            if string_2 not in fle:
                files.append(fle)
    return files


def locate_leaders(filename):
    donor = (re.split('%20|%3a|%29', filename))[8]
    sample = (re.split('%20|%3a|%29', filename))[10]
    output = (output_path + '%s_%s_locations.tsv') % (sample, donor)
    #query_list
    query_file = ('fishers_output_%s_%s_filteredfdr.txt') % (sample, donor)
    query_file = query_file.replace('onor', '')
    print query_file
    print output
    if query_file in unfinished:
        query_file = output_path + query_file
        query_list = []
        with open(query_file, 'r') as inF:
            next(inF)
            for line in inF:
                query = (line.split('\t'))[0]
                query_list.append(query)
        filename = path + filename
        bamfiles_list = []
        bamfiles_list2 = []
        with open(filename, 'r') as inF:
            for line in inF:
                if 'VHE' in line:
                    line = line.split('\t')
                    if line[1] == '0':
                        bamfiles_list.append(line[9])
                        bamfiles_list2.append([line[9], line[2], (str(line[3]) + ',+')])
                    elif line[1] == '16':
                        bamfiles_list.append(line[9])
                        bamfiles_list2.append([line[9], line[2], (str(line[3]) + ',-')])
                    else:
                        pass
        o = open(output, 'ab+')
        for i in range(0, len(query_list)):
            value = list(map(lambda x: str(x).startswith(query_list[i]), (bamfiles_list)))
            indexes = ([y for y, j in zip(count(), value) if j == True])
            location = list(map(lambda x : str((bamfiles_list2[x])[1]) + ' ' + str((bamfiles_list2[x])[2]), indexes))
            total = str(len(location))
            o.write(str(query_list[i]) + '\t' + total)
            location = Counter(location)
            location = Counter.most_common(location)
            for l in location:
                o.write('\t' + str(l[0]) + ':' + str(l[1]))
            o.write('\n')
        return output


file_list = get_files(path, '.sam', 'mock')

# print file_list

#for input_file in file_list:
#    locate_leaders(input_file)

# part 1 get chromosomal locations
# this part take s along time and the time differs between samfile.
pool = mp.Pool(processes=16)
results = pool.map(locate_leaders, file_list)
