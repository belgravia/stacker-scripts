import sys, csv  # unfinished script sorry`

try:
	gtf = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py gtf outfilename\n')
	sys.exit(1)

lengths = {}  # isoform: length
for line in gtf:
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand, info = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8]
	if ty != 'exon':
		continue
	if strand == '+':
		if chrom not in annotpos:
			annotpos[chrom] = {}
		annotpos[chrom][start] = {}  # the start of an exon on the positive strand is a junction 3'
	elif strand == '-':
		if chrom not in annotmin:
			annotmin[chrom] = {}
		annotmin[chrom][end] = {}


with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
