import sys, csv

try:
    psl = open(sys.argv[1])
    mode = sys.argv[2]
    outfilename = sys.argv[3]
except:
    sys.stderr.write('usage: script.py pslfile [naive/nomm/identity] outfilename\n')
    sys.exit(1)

reads = {}
for line in psl:
    line = line.rstrip().split('\t')
    matches, transcript, tlen, read = int(line[0]), line[13], int(line[14]), line[9]
    if read not in reads:
		reads[read] = []
    if mode == 'naive':  # use with primary alignments only
    	reads[read] += [(transcript, 1)]
    elif mode == 'nomm':  # no multiple mappers. should supply psl with all sub alignments too
    	reads[read] += [(transcript, 1)]
    elif mode == 'identity':
	    reads[read] += [(transcript, matches/float(tlen))]
    else:
		sys.stderr.write('unrecognized mode\n')
		sys.exit()

quant = {}
for read in reads:
	if mode == 'nomm' and len(reads[read]) == 1:
		transcript = reads[read][0][0]
		if transcript not in quant:
			quant[transcript] = 0
		quant[transcript] += 1
	elif mode == 'naive':
		tot = len(reads[read])
		for t in reads[read]:
			if t[0] not in quant:
				quant[t[0]] = 0
			quant[t[0]] += t[1]/float(tot)
		# print(read, reads[read])
	elif mode == 'identity':
		totsum = sum([t[1] for t in reads[read]])
		for t in reads[read]:
			if t[0] not in quant:
				quant[t[0]] = 0
			quant[t[0]] += t[1]/float(totsum)

with open(outfilename, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    for transcript in quant:
    	writer.writerow([transcript, quant[transcript]])