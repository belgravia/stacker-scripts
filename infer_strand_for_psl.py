import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	print('usage: script.py psl gtf_v26 outfilename')
	sys.exit(1)

annotmin = {}
annotpos = {}
dups = 0
geneset = set()
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
	geneset.add(gene)
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

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	dist = {}
	prevline, maxdist = '', 50
	for line in psl:
		line = line.rstrip().split('\t')
		chrom  = line[13]
		starts = [int(x) for x in line[20].split(',')[:-1]]  # block/exon starts
		ends = [int(x) + y for x,y in zip(line[18].split(',')[:-1], starts)]

		jstarts = ends[:-1]  # junction/intron starts
		jends = starts[1:]

		if len(ends) == 1:
			jstarts = starts
			jends = ends


		strands = []  # tally of predicted strands for each junction
		gene_candidates = []
		alljunctionstrands = []
		for start, end in zip(jstarts, jends):
			if chrom not in annotpos or chrom not in annotmin:
				continue
			wiggle = find_wiggle(start, annotpos[chrom], annotmin[chrom], maxdist)
			if wiggle == maxdist:
				startstrand = '.'
				gene = ''
			elif start+wiggle in annotpos[chrom]:
				startstrand = '+'
				gene = annotpos[chrom][start+wiggle]
			else:
				startstrand = '-'
				gene = annotmin[chrom][start+wiggle]
			# if wiggle in dist:
			# 	dist[wiggle] += int(line[4])
			# else:
			# 	dist[wiggle] = int(line[4])
			wiggle = find_wiggle(end, annotpos[chrom], annotmin[chrom], maxdist)
			if wiggle == maxdist:
				endstrand = '.'
			elif end+wiggle in annotpos[chrom]:
				endstrand = '+'
			else:
				endstrand = '-'

			consensus = startstrand if startstrand == endstrand else '.'
			strands += [consensus]
			alljunctionstrands += [startstrand, endstrand]
			gene_candidates += [gene]

		if strands.count('+') > strands.count('-'):
			consensus = '+'
			gene = gene_candidates[strands.index('+')]
		elif strands.count('-') > strands.count('+'):
			consensus = '-'
			gene = gene_candidates[strands.index('-')]
		else:
			if alljunctionstrands.count('+') > alljunctionstrands.count('-'):
				consensus = '+'
			elif alljunctionstrands.count('-') < alljunctionstrands.count('+'):
				consensus = '-'
			else:
				consensus = '.'
				gene = chrom+':'+str(jstarts[0])




		line[8] = consensus
		line[9] += '_' + gene
		writer.writerow(line)
# print(dist)
# print(dups)