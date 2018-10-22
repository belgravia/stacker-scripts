import sys, csv

try:
    tsv = open(sys.argv[1])
    outfilename = sys.argv[2]
except:
    print('usage: script.py tsv outfilename')
    sys.exit()

counts = {}
tsv.readline()  # header
for line in tsv:
    line = line.rstrip().split('\t')
    gene = line[0][line[0].find('ENSG'):]
    gene = gene[:gene.find('|')]
    if gene not in counts:
        counts[gene] = 0
    counts[gene] += float(line[3])

with open(outfilename, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    for g in counts:
        writer.writerow([g, counts[g]])
