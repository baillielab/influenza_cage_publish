infile = 'leader_counts.txt'

leader = []

with open(infile, 'r') as inF:
	for line in inF:
		linea = line.split('\t')
		leader.append(linea[0])

print leader[0:10]

print len(leader)
print len(list(set(leader)))