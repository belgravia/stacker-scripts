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
bad=False
dups = 0
geneset = set()
for line in gtf:   # i really should write a gtf parser..
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand, gene = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8].split(';')[4]
	if ty != 'exon':
		continue
	gene = gene[len('gene_name') + 3:-1]
	geneset.add(gene)
	if strand == '+':
		if chrom not in annotpos:
			annotpos[chrom] = {}
		if start in annotpos[chrom] and gene != annotpos[chrom][start]:
			dups += 1
			# annotpos[chrom][start] += [gene]
		else:
			annotpos[chrom][start] = gene
		if end in annotpos[chrom] and gene != annotpos[chrom][end]:
			dups += 1
			# annotpos[chrom][end] += [gene]
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

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	dist = {}
	prevline = ''
	for line in psl:
		line = line.rstrip().split('\t')
		chrom  = line[13]
		starts = [int(x) for x in line[20].split(',')[:-1]]
		ends = [int(x) + y for x,y in zip(line[18].split(',')[:-1], starts)]

		jstarts = ends[:-1]
		jends = starts[1:]

		print(chrom, starts, ends)
		strands = []
		gene_candidates = []
		for start, end in zip(jstarts, jends):
			wiggle = 0
			if chrom not in annotpos or chrom not in annotmin:
				continue
			while start + wiggle not in annotpos[chrom] and start + wiggle not in annotmin[chrom]:
				if wiggle == 100:
					break
				if wiggle == 0:
					wiggle += 1
				elif wiggle >= 0:
					wiggle = wiggle * -1
				else:
					wiggle = (wiggle-1) * -1
			if wiggle == 100:
				startstrand = '.'
				gene = ''
			elif start+wiggle in annotpos[chrom]:
				startstrand = '+'
				gene = annotpos[chrom][start+wiggle]
			else:
				startstrand = '-'
				gene = annotmin[chrom][start+wiggle]
			if wiggle in dist:
				dist[wiggle] += int(line[4])
			else:
				dist[wiggle] = int(line[4])
			wiggle = 0
			while end + wiggle not in annotpos[chrom] and end + wiggle not in annotmin[chrom]:
				if wiggle == 100:
					break
				if wiggle == 0:
					wiggle += 1
				elif wiggle >= 0:
					wiggle = wiggle * -1
				else:
					wiggle = (wiggle-1) * -1
			if wiggle == 100:
				endstrand = '.'
			elif end+wiggle in annotpos[chrom]:
				endstrand = '+'
			else:
				endstrand = '-'

			consensus = startstrand if startstrand == endstrand else '.'
			strands += [consensus]
			gene_candidates += [gene]

		if gene_candidates:
			print(gene_candidates)
		print('strands', strands)
		if strands.count('+') > strands.count('-'):
			consensus = '+'
			gene = gene_candidates[strands.index('+')]
		elif strands.count('-') > strands.count('+'):
			consensus = '-'
			gene = gene_candidates[strands.index('-')]
		else:
			consensus = '.'
			gene = ''
		line[8] = consensus
		line[9] += '_' + gene
		# if prevline == '':
			# prevline = line
			# continue
		# elif line[:3] != prevline[:3]:
			# writer.writerow(prevline)
			# prevline = line
		# else:
			# prevline[4] = int(line[4]) + int(prevline[4])
		writer.writerow(line)
print(dist)
print(dups)
if bad:
	print('bad')