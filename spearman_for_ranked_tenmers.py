import scipy.stats
import os
from scipy.stats import spearmanr
import sys

def Spearman_segments(rankedlist, dict1, dict2, name_file_1, name_file_2, output):
    o = open(outfile, 'ab')
    len_1 = len(dict1)
    len_2 = len(dict2)
    new_1 = []
    new_2 = []
    i = 0
    while i < len(rankedlist):
        leader = rankedlist[i]
        if leader in dict1:
            rank = dict1[leader]
            new_1.append(rank)
        else:
            new_1.append(len_1)

        if leader in dict2:
            rank = dict2[leader]
            new_2.append(rank)
        else:
            new_2.append(len_2)
        i = i + 1
    rho, p = scipy.stats.spearmanr(new_1,new_2)
    outlist = [name_file_1, name_file_2, str(rho), str(p), '\n']
    output = '\t'.join(outlist)
    o.write(output)
    o.close()

def get_ranked_leaders(ranked):
    ranked_leaders = []
    with open(ranked, 'r') as inF:
        for line in inF:
            ranked_leaders.append((line.split('\t'))[0])
    return ranked_leaders


def get_for_intersect(file1):
    lst1 = []
    with open(file1, 'r') as inF:
        for line in inF:
            line = line.split('\t')
            ten = line[0]
            lst1.append(ten)
    return lst1


def Intersection(lst1, lst2): 
    '''Python program to illustrate the intersection. 
       of two lists using set() and intersection()
       Rerurns overlap of lists.
    '''
    out = list((set(lst1) & set(lst2)))
    return out
      
#lst1 = ['dog', 'cat', 'possum', 'fish', 'pig', 'degu']
#lst2 = ['dog', 'cat', 'degu', 'fish', 'cow', 'goat']
#print(Intersection(lst1, lst2)) 

def create_dictionaries_for_spearmans(overlapList, file1, file2):
    dict1 = {}
    with open(file1) as inF:
        i = 1
        for line in inF:
            line = line.split('\t')
            if line[0] in overlapList:
                dict1[line[0]] = i
                i = i + 1

    dict2 = {}
    with open(file2) as inF:
        next(inF)

        tenmer_list = []
        for line in inF:
            line = line.split('\t')
            ten = line[0]
            perc = float(line[1])
            if ten in overlapList:
                tenmer_list.append([ten, perc])
        output = sorted(tenmer_list, reverse = True, key=lambda x: float(x[1]))
        print output[:10]
        i = 1
        while i < len(output):
            dict2[(output[i])[0]] = i
            i = i + 1
    name_file_1 = (file1.split('.'))[0]
    name_file_2 = (file2.split('.'))[0]
    return dict1, dict2, name_file_1, name_file_2

def find_files(path, string):
    import os
    files = []
    dirs = os.listdir(path)
    for fle in dirs:
        if string in fle:
            files.append(fle)
    return files

rankedlist = get_ranked_leaders('ranked_tenmers_all_samples.txt')
current = os.getcwd()
#
#file1 = 'ranked_tenmers_all_samples.txt'

f = sys.argv[1]
outfile = (f.split('.'))[0] + '_spearman.txt'
lst2 = get_for_intersect(f)
files = find_files(current, 'Segment')

for file1 in files:
    if file1 == f:
        pass
    else:
        lst1 = get_for_intersect(file1)
        overlapList = Intersection(lst1, lst2) #finds only overlapping elements of the lists
        dict1, dict2, name_file_1, name_file_2 = create_dictionaries_for_spearmans(overlapList, file1, f)
        Spearman_segments(rankedlist, dict1, dict2, name_file_1, name_file_2, outfile)
