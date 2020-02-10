import random

def ENTdictionary(sourcefile):
    ent_genename_dictionary = {}
    with open(sourcefile, 'r') as inF:
        next(inF)
        for line in inF:
            linea = line.split('\t')
            ent_genename_dictionary[linea[1]] = linea[2].strip()

    return ent_genename_dictionary

gene_names = ENTdictionary('grch37_180401_geneid_transid_genename_mart_export.txt')

infile = 'promoter_list_unfiltered.tsv'
i = 6
while i <= 100:
    #create promoter list file
    check_list = []
    #remove_output
    outfile = (('promoter_list_filtered_%s.tsv')%(str(i)))
    o = open(outfile, 'w')
    outlist = ['tenmer', 'hour', 'donor', 
               'FDR', 'OR', 'snatched', 'unsnatched', 
               'promoter', 'single_prom', 'gene', '\n']
    output = '\t'.join(outlist)
    o.write(output)
    with open(infile, 'r') as inF:
        next(inF)
        for line in inF:
            linea = (line.strip()).split('\t')
            if len(linea) >= 8:
                if '..' in linea[7]:
                    outprom = linea[7]
                    outgene = (outprom.split('@'))[1]
                    outlist = linea[0:7] + [outprom, outgene, '\n']
                else:

                    promoter = linea[7].split(',')
                    if len(promoter) > 1:
                        single_prom = random.choice(promoter) #randomly choose a promoter
                    else:
                        single_prom = promoter[0]
                    gene = (single_prom.split('@'))[1]

                    if 'ENST' in gene:
                        if gene in gene_names:
                            gene = gene_names[gene]

                    outprom = (single_prom.split('@'))[0] + '@' + gene
                    outgene = gene
                    outlist = linea[0:7] + [outprom, outgene, '\n']
            
            else:
                outlist = linea[0:7] + ['NA', 'NA', '\n']
            
            output = '\t'.join(outlist)
            o.write(output)
    i = i + 1
