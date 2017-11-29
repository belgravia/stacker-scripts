import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	# outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py psl gtf\n')
	sys.exit(1)

prev_transcript = ''
annotated = {}
for line in gtf:  # extract all exons from the gtf, keep exons grouped by transcript
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
	if ty != 'exon':
		continue
	this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('gene_type')-3]  # p specific to gencode v24
	# this_transcript = line[8][line[8].find('transcript_id')+15:line[8].find('exon_assignment')-3]  # p specific to sirvs
	if this_transcript != prev_transcript:
		if prev_transcript:
			if transcript_chrom not in annotated:
				annotated[transcript_chrom] = {}
			annotated[transcript_chrom][prev_transcript] = {}
			annotated[transcript_chrom][prev_transcript]['sizes'] = blocksizes
			annotated[transcript_chrom][prev_transcript]['starts'] = blockstarts
		prev_transcript = this_transcript
		transcript_start = start  # convert to 0-index
		transcript_chrom = chrom 
		blockstarts = [int(start) - 1]
		blocksizes = [int(end) - int(start) + 1]
	else:
		blockstarts += [int(start) - 1]
		blocksizes += [int(end) - int(start) + 1]

novel = 0
total = 0
for line in psl:
	line = line.rstrip().split('\t')
	chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])
	blocksizes = [int(x) for x in line[18].split(',')[:-1]]
	blockstarts = [int(x) for x in line[20].split(',')[:-1]]
	found = False
	total += 1
	if chrom not in annotated:
		continue
	for transcript in annotated[chrom]:
		astarts = annotated[chrom][transcript]['starts']  # annotated starts
		asizes = annotated[chrom][transcript]['sizes']
		if asizes == blocksizes and astarts == blockstarts or \
			(len(blocksizes) > 1 and len(astarts) == len(blockstarts) and \
			blockstarts[1:] == astarts[1:] and blocksizes[1:-1] == asizes[1:-1] and \
			blockstarts[0] + blocksizes[0] == astarts[0] + asizes[0]): # classify isoforms that only differ by tes and tss as the same
		   found = True
		   sys.stderr.write('Same isoform: {}, {}\n'.format(transcript, name))
		   novel += 1
		   break
sys.stderr.write('{} out of {} isoforms are novel\n'.format(novel, total))

