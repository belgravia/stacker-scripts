import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	# outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py psl gtf\n')
	sys.exit(1)

def get_junctions(line):
	junctions = set()
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) + 1 for n in line[18].split(',')[:-1]]  # add 1 for indexing pupropses
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		# junctions += [(starts[b]+sizes[b], starts[b+1])]
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

prev_transcript = ''
annotated_juncs = {}
# for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
# 	if line.startswith('#'):
# 		continue
# 	line = line.rstrip().split('\t')
# 	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
# 	if ty != 'exon':
# 		continue
# 	this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('gene_type')-3]  # p specific to gencode v24
# 	# this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('exon_assignment')-3]  # p specific to sirvs
# 	if this_transcript != prev_transcript:
# 		if prev_transcript:
# 			if transcript_chrom not in annotated:
# 				annotated[transcript_chrom] = {}
# 			annotated[transcript_chrom][prev_transcript] = {}
# 			annotated[transcript_chrom][prev_transcript]['sizes'] = blocksizes
# 			annotated[transcript_chrom][prev_transcript]['starts'] = blockstarts
# 		prev_transcript = this_transcript
# 		transcript_start = start  # convert to 0-index
# 		transcript_chrom = chrom 
# 		blockstarts = [int(start) - 1]
# 		blocksizes = [int(end) - int(start) + 1]
# 	else:
# 		blockstarts += [int(start) - 1]
# 		blocksizes += [int(end) - int(start) + 1]

for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'exon':
		continue
	if chrom not in annotated_juncs:
		annotated_juncs[chrom] = []
	this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('gene_type')-3]  # p specific to gencode v24
	this_transcript = line[8][line[8].find('transcript_id')+15:]
	this_transcript = this_transcript[:this_transcript.find('"')]

	# this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('exon_assignment')-3]  # p specific to sirvs
	if this_transcript != prev_transcript:
		if prev_transcript:
			annotated_juncs[chrom] += [junctions]
			# print(junctions)
			# print(prev_transcript)
			# break
		junctions = set()
		prev_transcript = this_transcript
	else:
		junctions.add((prev_end, start))
	prev_end = end

novel, total = 0, 0
for line in psl:
	line = line.rstrip().split('\t')
	chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])
	junctions = get_junctions(line)
	# print(chrom)
	print(junctions)
	sys.exit()
	total += 1
	if junctions not in annotated_juncs[chrom]:
		print('\t'.join([str(l) for l in line]))
		novel += 1
sys.stderr.write('{} out of {} isoforms are novel\n'.format(novel, total))

