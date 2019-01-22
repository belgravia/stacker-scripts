import sys, csv, re
import pysam

try:
	if sys.argv[1][-3:] == 'bam':
		sam = pysam.AlignmentFile(sys.argv[1], 'rb')
	else:
		sam = pysam.AlignmentFile(sys.argv[1], 'r')
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py sam/bam outtxt\n')
	sys.stderr.write('outputs junction (or no junction) for every read\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in sam:
		if not line.cigarstring or line.cigarstring == '*':
			writer.writerow([line.query_name, 1])
			continue
		intron = sam.find_introns(l for l in [line]) # test
		if not intron:
			intron = 0
		else:
			intron = 'chr' + line.reference_name + ':' + '-'.join([str(x) for x in intron.keys()[0]])
		writer.writerow([line.query_name, intron])
		continue


		if line[0] == '@':
			continue
		# line = line.rstrip().split('\t')
		qname, flag, tname, pos, cigar = line[0], int(line[1]), line[2], int(line[3]), line[5]
		if tname == '*':
			continue
		pos = pos - 1
		matches = re.findall('([0-9]+)([A-Z])', cigar)
		matchlen = mismatches = relstart = qstart = qconsumed = 0
		blocksizes, relblockstarts, qstarts = [], [], []
		tend = pos
		qnuminsert = 0
		qbaseinsert = 0
		tnuminsert = 0
		tbaseinsert = 0  # deletion
		qsize_backup = 0
		junction = ''
		for m in matches:
			num, op = int(m[0]), m[1]
			if op == 'M':  # consumes reference
				blocksizes += [num]
				relblockstarts += [relstart]
				qstarts += [qconsumed]
				relstart += num
				matchlen += num
				tend += num
				qconsumed += num
				qsize_backup += num
			elif op == 'D':  # consumes reference
				relstart += num
				mismatches += num
				tend += num
				qnuminsert += num
			elif op == 'I':
				qconsumed += num
				tbaseinsert += num
				tnuminsert += 1
				qsize_backup += num
			elif op == 'N':  # consumes reference
				tend += num
				junction = tname +':'+ str(relstart+pos) +'-'+ str(relstart+num+pos)
				break
				relstart += num
			elif op == 'S':
				if not qstart and not matchlen:
					qstart = num
				qsize_backup += num
			elif op == 'H':
				if not qstart and not matchlen:
					qstart = num
					relstart += num
				qsize_backup += num  # technically does not consume q but useful when comparing a read's secondary alignments
			else:
					sys.stderr.write(op + '\n')
		qend = qconsumed + qstart
		ncount = seq.count('N')
		qsize = len(seq)
		qsize = qsize_backup
		tsize = chromsizes[tname]  # chromosome length
		tstart = pos
		strand = '-' if flag & 0x10 else '+'  # flag&0x10 is 1 when the strand is -
		blockstarts = [str(pos + s) for s in relblockstarts]
		blockcount = len(blockstarts)
		qstarts = ','.join([str(qstart + s) for s in qstarts]) + ','
		blocksizes = ','.join([str(s) for s in blocksizes]) + ','
		relblockstarts = ','.join([str(s) for s in relblockstarts]) + ','
		blockstarts = ','.join(blockstarts) + ','
		writer.writerow([qname, junction])