import sys, csv

try:
	txt = open(sys.argv[1])
	fasta = open(sys.argv[2])
	outfilename5 = sys.argv[3]
	outfilename3 = sys.argv[4]
	if len(sys.argv)>5:
		window = int(sys.argv[5])
	else:
		window = 10
	if len(sys.argv)>6:
		gtf = open(sys.argv[6])
		gtfused = True
		maxdist = 100
	else:
		gtf = ''
		gtfused = False
except:
	sys.stderr.write('usage: script.py coord.txt genome.fasta out5.fa out3.fa [num_bp_around] [gtf]\n')
	sys.stderr.write('txtfile format should be: chr start end strand\n')
	sys.exit(1)

def revcomp(seq):
	seq = seq.upper()
	seq = seq.replace('A', 'X').replace('T', 'A').replace('G', 'Y').replace('C', 'G')
	seq = seq.replace('X', 'T').replace('Y', 'C')
	return seq[::-1]

annotmin = {}
annotpos = {}
dups = 0
if gtf:
	for line in gtf:
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end, strand, gene = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8]
		if ty != 'exon':
			firstexon = True
			continue
		gene = gene[gene.find('gene_id')+len('gene_id')+2:]
		gene = gene[:gene.find('"')]
		if strand == '+':
			if chrom not in annotpos:
				annotpos[chrom] = {}
			if start in annotpos[chrom] and gene != annotpos[chrom][start]:
				dups += 1
			else:
				annotpos[chrom][start] = gene
			if end in annotpos[chrom] and gene != annotpos[chrom][end]:
				dups += 1
			else:
				annotpos[chrom][end] = gene
		else:
			if chrom not in annotmin:
				annotmin[chrom] = {}
			if start in annotmin[chrom] and gene != annotmin[chrom][start]:
				dups += 1
			else:
				annotmin[chrom][start] = gene
			if end in annotmin[chrom] and gene != annotmin[chrom][end]:
				dups += 1
			else:
				annotmin[chrom][end] = gene
		firstexon = False

def find_wiggle(coord, annot, annot2={}, maxdist=100):
	""" Finds the distance between coordinate and the closest annotated pos in annot dict. """
	wiggle = 0
	while coord + wiggle not in annot and coord + wiggle not in annot2:
		if wiggle == maxdist:
			break
		if wiggle == 0:
			wiggle += 1
		elif wiggle >= 0:
			wiggle = wiggle * -1
		else:
			wiggle = (wiggle-1) * -1
	return wiggle

query = {}
for line in txt:
	line = line.rstrip().split('\t')
	chrom, pos1, pos2, strand = line[:4]
	# chrom, pos1, pos2, pval, strand = line[:5]
	# pos2 = line[11]  # just make sure it's grabbing the right columns in the txtfile
	# if float(pval) > 0.1:
	# 	continue
	pos1, pos2 = int(pos1), int(pos2)
	if chrom not in query:
		query[chrom] = {}
	if strand not in query[chrom]:
		query[chrom][strand] = {}
	if 'left' not in query[chrom][strand]:
		query[chrom][strand]['left'] = []
		query[chrom][strand]['right'] = []
	if int(pos2) < int(pos1):
		pos1, pos2 = pos2, pos1
	if gtfused:
		wiggle = find_wiggle(pos1, annotpos[chrom], annotmin[chrom], maxdist)
		if wiggle != maxdist:
			pos1 += wiggle
		wiggle = find_wiggle(pos2, annotpos[chrom], annotmin[chrom], maxdist)
		if wiggle != maxdist:
			pos2 += wiggle
	query[chrom][strand]['left'] += [pos1]
	query[chrom][strand]['right'] += [pos2]

for chrom in query:
	for strand in query[chrom]:
		for lr in query[chrom][strand]:
			query[chrom][strand][lr] = sorted(query[chrom][strand][lr])

with open(outfilename3, 'wt') as outfile3, open(outfilename5, 'wt') as outfile5:
	writer5 = csv.writer(outfile5, delimiter='\t')
	writer3 = csv.writer(outfile3, delimiter='\t')
	seq, chrom = '', ''
	for line in fasta:
		line = line.rstrip()
		if line.startswith('>'):
			if not chrom:
				chrom = line[1:]
				continue
			if chrom in query:
				for strand in query[chrom]:
					for lr in query[chrom][strand]:
						for pos in query[chrom][strand][lr]:
							if strand == '+' and lr == 'left':
								writer5.writerow([seq[pos-window:pos+window].upper()])  # formerly pos:pos+2
							elif strand == '+' and lr == 'right':
								writer3.writerow([seq[pos-window:pos+window].upper()])
							elif strand == '-' and lr == 'left':
								writer3.writerow([revcomp(seq[pos-window:pos+window])])
							else:
								writer5.writerow([revcomp(seq[pos-window:pos+window])])
			chrom = line[1:]
		seq = line.rstrip()

