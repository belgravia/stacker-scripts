import sys, csv

try:
	psl = open(sys.argv[1])
	fastq = open(sys.argv[2])
	if len(sys.argv) == 4:
		ishoulduseargparse = 4/0
	if len(sys.argv) > 4:
		mutchrom = sys.argv[3]
		gpos = int(sys.argv[4])
	else:
		mutchrom = 'chr2'
		gpos = 197402110  # position of mutation of interest
except:
	sys.stderr.write('usage: script.py pslfile fastq [chr] [position]\n')
	sys.stderr.write('assumes that the gene is on the minus strand\n')
	sys.exit(1)


aa_dict = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"Z", "UAG":"Z",  # using Z for stop
    "UGU":"C", "UGC":"C", "UGA":"Z", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",} 
def revcomp(seq):
	seq = seq.replace('A', 'X').replace('T', 'A').replace('X', 'T')
	seq = seq.replace('G', 'X').replace('C', 'G').replace('X', 'C')
	return seq[::-1]
def translate_seq(seq):
	aa_string = ''
	for i in range(0,len(seq),3):
		if i + 3 > len(seq):
			break
		aa_string += aa_dict[seq[i:i+3]]
	return aa_string

read_info = {}
for line in psl:
	line = line.rstrip().split('\t')
	strand, read, chrom = line[8], line[9][:-4], line[13]
	if chrom != mutchrom:
		continue
	blocksizes = [int(x) for x in line[18].split(',')[:-1]]
	blockstarts = [int(x) for x in line[20].split(',')[:-1]]
	if strand == '+':
		blockstarts_query = [int(x) for x in line[19].split(',')[:-1]]
		for b in range(len(blockstarts)):
			if blockstarts[b] <= gpos and blockstarts[b] + blocksizes[b] > gpos:
				read_info[read] = (blockstarts_query[b] + gpos - blockstarts[b] - 3, 
					strand)
	else:
		blockstarts = blockstarts[::-1]
		blocksizes = blocksizes[::-1]
		querypos = 0
		for b in range(len(blockstarts)):
			if blockstarts[b] < gpos and blockstarts[b] + blocksizes[b] >= gpos:
				read_info[read] = (querypos + blockstarts[b] + blocksizes[b] - gpos + 3,
					strand)
			querypos += blocksizes[b]

num = 0
lookup = False
for line in fastq:
	if num % 4 == 0:
		read = line.split()[0][1:-4]
		if read in read_info:
			lookup = True
	elif lookup and num % 4 == 1:
		pos, strand = read_info[read]
		codon = line[pos:pos+9]
		if strand == '-':  # change these lines depending on what strand the gene is on
			longercodon = line[pos:pos+18]
			longercodon = line.rstrip()
		else:
			codon = revcomp(codon)
			longercodon = revcomp(line[pos-6:pos+12])
		rnacodon = longercodon.replace('T', 'U')
		print('{} {} {}\t{}\t{}\t{}\n'.format(read, codon, pos+9, translate_seq(rnacodon), 
			translate_seq(rnacodon[1:]), translate_seq(rnacodon[2:]))) 
		if 'b3f0' in read or '83397' in read or '489d' in read:
			print('{} {}'.format(read, rnacodon))
		lookup = False
	num += 1
