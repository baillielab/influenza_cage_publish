

def create_dict(tenmer):
    loc_dict = {}
    with open(direct + 'location_table_full.5.10.CR.tsv', 'r') as inF:
        for line in inF:
            if tenmer in line:
                linea = line.split('|')
                chrom = linea[1]
                start = linea[2]
                prom = linea[6]
                loc_dict[chrom+'_'+start] = prom
    return loc_dict


thisseq_snatched = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/source_data/10mercount_tables/thisseq_snatched'
tenmers = []
with open(thisseq_snatched, 'r') as inF:
    for line in inF:
        if 'seq' in line:
            pass
        else:
            linea = line.split('\t')
            tenmers.append(linea[0])
collated_tenmers = 'collated_tenmers_locations.tsv'

f = open(collated_tenmers, 'r')
lines = f.readlines()
f.close()

direct = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/4_leaders/leader_tables/2_name/'
outfile = 'promoter_list.tsv'
o = open(outfile, 'w')

for tenmer in tenmers:
    loc_dict = create_dict(tenmer)
    tenmer_lines = []
    for line in lines:
        if tenmer in line:
            tenmer_lines.append(line)
    for entry in tenmer_lines:
        entry = entry.split('\t')
        locations = (entry[-2]).split(',')
        loc_list = []
        for loc in locations[1:]:
            if 'chr' in loc:
                chrom = loc.replace("'", "")
                chrom = chrom.replace('"', '')
                chrom = chrom.replace('[', '')
                chrom = chrom.replace(' c', 'c')
                chrom = chrom.replace(' ', '_')
                loc_list.append(chrom)
        prom_list = []
        for loc in loc_list:
            try:
                promoter = loc_dict[loc]
                prom_list.append(promoter)
            except KeyError:
                pass
        if len(prom_list) > 0:
            prom_list = ','.join(prom_list)
            outline = [entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], entry[6], str(prom_list), '\n']
        else:
            outline = [entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], entry[6], ' ', '\n']
 
        output = '\t'.join(outline)
        o.write(output)