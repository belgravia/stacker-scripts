import sys, csv

try:
	gtf = open(sys.argv[1])
	genomefa = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py gtf genome_fa dists_outfile\n')
	sys.exit()

def first_non_gag_trimer(seq, strand, query='CT', bump=0):
	""" Strand is essentially the search direction. -:L->R, +:R->L """
	if len(seq) < 3:
		return -1
	if strand == '-':
		pos = seq.find(query)
		if pos < 0:  # not found
			return pos
		if len(seq) <= pos+2:
			return -1
		if seq[pos+2] != 'C':
			return 101 - len(seq[pos:]) + 1
		return first_non_gag_trimer(seq[pos+2:], strand)  # keep searching until non CTC
	else:
		pos = seq.rfind(query)
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
				pos = first_non_gag_trimer(seq[tp+2:tp+103], '-')  # 100 bp into intron
				pos_other = first_non_gag_trimer(seq[tp-101:tp],'+') # 100 bp into exon
				print(seq[tp+2:tp+103], pos, strand)
				print(seq[tp-101:tp], pos_other, strand)

				if pos_other != -1 and pos_other < pos:
					pos = pos_other
					toneg = True
				else:
					toneg = False
				if pos >= 0:
					if toneg:
						pos = -1 * pos
					print(pos)
					writer.writerow([chrom, tp, pos, '-'])

			for tp in ss[chrom]['+']:
				pos = first_non_gag_trimer(seq[tp-101:tp], '+', 'AG')
				pos_other = first_non_gag_trimer(seq[tp+2:tp+103], '-', 'AG')
				# print(seq[tp-101:tp], pos, '+')
				if pos_other - tp < tp - pos:
					pos = pos_other
				if pos >= 0:
					writer.writerow([chrom, tp, pos, '+'])
			chrom = line[1:]
			seq = ''
			continue
		else:
			seq += line
