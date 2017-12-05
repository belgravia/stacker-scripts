import sys, csv, re

try:
	sam = open(sys.argv[1])
	chromsizefile = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py samfile chromsizes outpsl\n')
	sys.stderr.write('written for minimap sams\n')
	sys.exit(1)

def cigar_to_blocks(matches):  # parses cigar string matches, returns columns 19 and 20 of psl
	return 

chromsizes = {}
for line in chromsizefile:
	line = line.rstrip().split('\t')
	chromsizes[line[0]] = line[1]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in sam:
		if line.startswith('@'):
			continue
		line = line.rstrip().split('\t')
		tname, flag, qname, pos, cigar, seq, qual = line[0], int(line[1]), line[2], int(line[3]), line[5], line[9], line[10]

		matches = re.findall('([0-9]+)([A-Z])', cigar)
		matchlen = mismatches = 0
		blocksizes, relblockstarts = [], []
		relstart = 0
		for m in matches:
			# if m[1] in ['M', 'D']:  # both matches and deletions consume reference
				# matchlen += int(m[0])
			num, op = int(m[0]), m[1]
			if op == 'M':
				blocksizes += [num]
				relblockstarts += [relstart]
				relstart += num
				matchlen += num
			elif op == 'D':
				relstart += num
				mismatches = 0


		ncount = seq.count('N')
		qsize = len(seq)
		qstart = 0
		qend = 0
		if qname == '*':
			continue
		tsize = chromsizes[qname]  # chromosome length
		tstart = 0
		tend = 0
		strand = '-' if flag & 0x10 else '+'  # flag&0x10 is 1 when the strand is -
		qnuminsert = 0
		qbaseinsert = 0
		tnuminsert = 0
		tbaseinsert = 0  # deletion
		blockstarts = [str(pos + s) for s in relblockstarts]
		blockcount = len(blockstarts)
		blocksizes = ','.join([str(s) for s in blocksizes]) + ','
		relblockstarts = ','.join([str(s) for s in relblockstarts]) + ','
		blockstarts = ','.join(blockstarts)
		writer.writerow([matchlen, mismatches, 0, ncount, qnuminsert, qbaseinsert, \
			tnuminsert, tbaseinsert, strand, qname, qsize, qstart, qend, \
			tname, tsize, tstart, tend, blockcount, blocksizes, relblockstarts, blockstarts])