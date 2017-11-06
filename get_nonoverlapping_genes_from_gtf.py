import sys, csv

try:
	gtf = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	print('usage: script.py gtf outfilename [smallgenesfile]')
	sys.exit(1)

ensg = set()
small = set()
gene_name = ''
numexons = 0
for line in gtf:   # i really should write a gtf parser..
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand, gene = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8].split(';')[0]
	if ty == 'gene':
		if numexons > 1 and gene_name: # excluding small RNAs, sorry. also things that aren't spliced sep10 2017
			ensg.add(gene_name)
		if end - start > 200:
			gene_name = gene[len('gene_id')+2:-1]
		else:
			gene_name = ''
		numexons = 0
	elif ty == 'transcript':
		numexons = 0
	elif ty == 'exon':
		numexons += 1

# sys.stderr.write('finished reading gtf\n')

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for gene in ensg:
		writer.writerow([gene])

if len(sys.argv) > 3:
	smallgenesfile = sys.argv[3]
	with open(smallgenesfile, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t')
		for gene in small:
			writer.writerow([gene])