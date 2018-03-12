import sys, csv

try:
	txt = open(sys.argv[1])
	fasta = open(sys.argv[2])
	outfilename5 = sys.argv[3]
	outfilename3 = sys.argv[4]
	if len(sys.argv)>5:
		window = int(sys.argv[5])
	else:
		window = 10
except:
	sys.stderr.write('usage: script.py coord.txt genome.fasta out5.txt out3.txt [num_bp_around]\n')
	sys.stderr.write('chr start end strand')
	sys.exit(1)

def revcomp(seq):
	seq = seq.upper()
	seq = seq.replace('A', 'X').replace('T', 'A').replace('G', 'Y').replace('C', 'G')
	seq = seq.replace('X', 'T').replace('Y', 'C')
	return seq[::-1]



query = {}
for line in txt:
	line = line.rstrip().split('\t')
	# chrom, pos1, pos2, pval, strand = line[:5]
	chrom, pos1, pos2, strand = line[:4]
	# pos2 = line[11]  # just make sure it's grabbing the right columns in the txtfile
	# if float(pval) > 0.05:
		# continue
	if chrom not in query:
		query[chrom] = {}
	if strand not in query[chrom]:
		query[chrom][strand] = {}
	if 'left' not in query[chrom][strand]:
		query[chrom][strand]['left'] = []
		query[chrom][strand]['right'] = []
	if int(pos2) < int(pos1):
		pos1, pos2 = pos2, pos1
	query[chrom][strand]['left'] += [int(pos1)]
	query[chrom][strand]['right'] += [int(pos2)]

for chrom in query:
	for strand in query[chrom]:
		for lr in query[chrom][strand]:
			query[chrom][strand][lr] = sorted(query[chrom][strand][lr])

with open(outfilename3, 'wt') as outfile3, open(outfilename5, 'wt') as outfile5:
	writer5 = csv.writer(outfile5, delimiter='\t')
	writer3 = csv.writer(outfile3, delimiter='\t')
	seq, chrom = '', ''
	for line in fasta:
		line = line.rstrip()
		if line.startswith('>'):
			if not chrom:
				chrom = line[1:]
				continue
			if chrom in query:
				for strand in query[chrom]:
					for lr in query[chrom][strand]:
						for pos in query[chrom][strand][lr]:
							if strand == '+' and lr == 'left':
								writer5.writerow([seq[pos-window:pos+window].upper()])  # formerly pos:pos+2
							elif strand == '+' and lr == 'right':
								writer3.writerow([seq[pos-window:pos+window].upper()])
							elif strand == '-' and lr == 'left':
								writer3.writerow([revcomp(seq[pos-window:pos+window])])
							else:
								writer5.writerow([revcomp(seq[pos-window:pos+window])])
			chrom = line[1:]
		seq = line.rstrip()

