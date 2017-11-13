import sys, csv

try:
	gtf = open(sys.argv[1])
	psl = open(sys.argv[2])
	if len(sys.argv) > 3:
		outfilename = sys.argv[3]
	else:
		outfilename = ''
except:
	sys.stderr.write('usage: script.py gtf psl [outfile]\n')
	sys.stderr.write('the out is for an annotated intron length distribution file\n')
	sys.exit(1)

def overlap(coords0, coords1, tolerance=0):  # takes two tuples of (start, end) coordinates
	if not tolerance:
		return coords1[0] >= coords0[0] and coords1[0] <= coords0[1] or \
			coords1[1] >= coords0[0] and coords1[1] <= coords0[1]  # partial overlap is sufficient
	return coords1[0]-tolerance <= coords0[0] and coords1[1]+tolerance >= coords0[1]  # complete coverage

annotated = {}
prev_transcript = ''
for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'exon':
		continue
	this_transcript = line[8][line[8].find('transcript_id')+2:line[8].find('gene_type')-3]  # p specific to gencode v24
	if this_transcript != prev_transcript:
		if prev_transcript:
			if transcript_chrom not in annotated:
				annotated[transcript_chrom] = {}
			annotated[transcript_chrom][prev_transcript] = {}
			annotated[transcript_chrom][prev_transcript]['sizes'] = blocksizes
			annotated[transcript_chrom][prev_transcript]['starts'] = blockstarts
			annotated[transcript_chrom][prev_transcript]['range'] = start, blockstarts[-1]+blocksizes[-1]
		prev_transcript = this_transcript
		transcript_start = start
		transcript_chrom = chrom
		blockstarts = [int(start)]
		blocksizes = [int(end) - int(start)]
	else:
		blockstarts += [int(start)]
		blocksizes += [int(end) - int(start)]

annotated_introns = {}
for chrom in annotated:  # extract all introns from the exon dictionary
	for transcript in annotated[chrom]:
		if chrom not in annotated_introns:
			annotated_introns[chrom] = set()
		sizes, starts = annotated[chrom][transcript]['sizes'], annotated[chrom][transcript]['starts']
		intron_start = starts[0]+sizes[0]
		for size, start in zip(sizes[1:], starts[1:]):
			annotated_introns[chrom].add((intron_start, start))
			intron_start = start + size

toremove = set()  # set of introns that have exons in another transcript that overlap
for chrom in annotated_introns:  # remove introns for which there exist exons, can only do so after reading in all exons
	# the only point of the annotated dictionary is so i don't have to read in the gtf twice
	for intron in annotated_introns[chrom]:
		for transcript in annotated[chrom]:
			sizes, starts = annotated[chrom][transcript]['sizes'], annotated[chrom][transcript]['starts']
			if not overlap(annotated[chrom][transcript]['range'], intron):
				continue
			for size, start in zip(sizes, starts):
				if overlap((start, start+size), intron[0], intron[1]):
					toremove.add(intron)
	annotated_introns[chrom] = annotated_introns[chrom] - toremove

isoforms = {}
for line in psl:
	line = line.rstrip().split('\t')
	chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])
	blocksizes = [int(x) for x in line[18].split(',')[:-1]]
	blockstarts = [int(x) for x in line[20].split(',')[:-1]]
	if chrom not in isoforms:
		isoforms[chrom] = {}
	isoforms[chrom][name] = {}
	isoforms[chrom][name]['entry'] = line
	isoforms[chrom][name]['sizes'] = blocksizes
	isoforms[chrom][name]['starts'] = blockstarts
	isoforms[chrom][name]['range'] = start, end
	isoforms[chrom][name]['ir'] = False  # detection of intron retention event

numunique = 0
intron_distribution = {}
for chrom in isoforms:
	for iname0 in isoforms[chrom]:
		for iname1 in annotated[chrom]:
			if iname0 == iname1:
				continue
			if not overlap(isoforms[chrom][iname0]['range'], annotated[chrom][iname1]['range']):
				continue
			starts0, sizes0 = isoforms[chrom][iname0]['starts'], isoforms[chrom][iname0]['sizes']
			starts1, sizes1 = annotated[chrom][iname1]['starts'], annotated[chrom][iname1]['sizes']
			prev5 = starts1[0]+sizes1[0]  # previous 5' end of isoform1's intron
			for start1, size1 in zip(starts1[1:], sizes1[1:]):  # change this part to be seeing if it's a unique intron
				if start1-prev5 not in intron_distribution:
					intron_distribution[start1-prev5] = 1
				else:
					intron_distribution[start1-prev5] += 1
				for start0, size0 in zip(starts0, sizes0):
					if start0 < prev5 + 10 and start0+size0 > start1-10:  # if isoform 0 has exon where isoform 1 has intron
						isoforms[chrom][iname0]['ir'] = True
						break
				if isoforms[chrom][iname0]['ir']:
					break
				prev5 = start1+size1


sys.stderr.write('# novel introns: {}\n'.format(numunique))

if not outfilename:
	sys.exit()

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	intron_sizes = sorted(intron_distribution.keys())
	for s in intron_sizes:
		writer.writerow([s, intron_distribution[s]])


