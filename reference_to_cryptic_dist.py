import sys, csv

try:
	gtf = open(sys.argv[1])
	genomefa = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py gtf genome_fa dists_outfile\n')
	sys.exit()

def first_non_gag_trimer(seq, strand, bump=0):
	if len(seq) < 3:
		return -1
	if strand == '-':
		pos = seq.find('CT')
		if pos < 0:  # not found
			return pos
		if len(seq) <= pos+2:
			return -1
		if seq[pos+2] != 'C':
			return 101 - len(seq[pos:]) + 1
		return first_non_gag_trimer(seq[pos+2:], strand)  # keep searching until non CTC
	else:
		pos = seq.rfind('AG')
		if pos < 0:
			return pos
		if pos == 0:
			return -1
		if seq[pos-1] != 'G':
			return 101 - pos - 2
		return first_non_gag_trimer(seq[:pos], strand)

ss = {}  # all annotated 3' SS
previd = ''
lastend = ''
for line in gtf:
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand, info = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8]
	transcript_id = info[info.find('transcript_id')+len('transcript_id')+2:]
	transcript_id = transcript_id[:transcript_id.find('"')]
	if chrom not in ss:
		ss[chrom] = {}
	if strand not in ss[chrom]:
		ss[chrom][strand] = set()
	if ty != 'exon':
		continue
	if transcript_id != previd or not previd:  # first exon
		lastend = end
		previd = transcript_id
		continue
	if strand == '-':
		ss[chrom][strand].add(lastend)# += [lastend]
	else:
		ss[chrom][strand].add(start)# += [start]
	lastend = end
	previd = transcript_id

# print(ss)
seq = ''
chrom = ''
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in genomefa:
		line = line.rstrip()
		if line.startswith('>'):
			if not seq:
				chrom = line[1:]
				continue
			if chrom not in ss:
				chrom = line[1:]
				seq = ''
				continue
			seq = seq.upper()
			for tp in ss[chrom]['-']: # threeprime
				pos = first_non_gag_trimer(seq[tp+2:tp+103], '-')
				# print(seq[tp:tp+101], pos, strand)
				if pos >= 0:
					print(seq[tp+2:tp+103])
					writer.writerow([chrom, tp, pos, '-'])
			for tp in ss[chrom]['+']:
				pos = first_non_gag_trimer(seq[tp-103:tp-2], '+')
				# print(seq[tp-101:tp], pos, '+')
				if pos >= 0:
					writer.writerow([chrom, tp, pos, '+'])
			chrom = line[1:]
			seq = ''
			continue
		else:
			seq += line
