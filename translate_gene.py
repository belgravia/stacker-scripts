import sys, csv

try:
	psl = open(sys.argv[1])
	genome = open(sys.argv[2])
	if len(sys.argv) == 4:
		ishoulduseargparse = 4/0
	if len(sys.argv) > 4:
		mutchrom = sys.argv[3]
		gpos = int(sys.argv[4])
	else:
		mutchrom = 'chr2'
		gpos = 197402110  # position of mutation of interest

except:
	sys.stderr.write('usage: script.py file.psl genome.fa\n')
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
		if 'N' in seq[i:i+3]:
			aa_string += 'O'
		else:
			aa_string += aa_dict[seq[i:i+3]]
	return aa_string[::-1]

def get_sequence(entry, seq):
	blocksizes = [int(x) for x in entry[18].split(',')[:-1]]
	blockstarts = [int(x) for x in entry[20].split(',')[:-1]]
	pulledseq = ''
	for block in range(len(blockstarts)):
		pulledseq += seq[blockstarts[block]:blockstarts[block]+blocksizes[block]]
	return pulledseq

psldata = {}
for line in psl:
	line = line.rstrip().split('\t')
	if line[9] in ['ENSG00000167996.15-5151_chr11_Isoform4_E_S_2', 
	'ENSG00000167996.15-5151_chr11_Isoform2_E_S_4', 'ENSG00000167996.15-5151_chr11_Isoform1_E_S_160']:

	# all the way across, 

	# if '8764f77e' in line[9] or 'baa3a9e2' in line[9] or \
	#   'fa84730e' in line[9]:  # third intron retained, all introns retained with small gap, 
	# # all introns retained
		print(line)
		if line[13] in psldata:
			psldata[line[13]] += [line]
		else:
			psldata[line[13]] = [line]

seq, chrom = '', ''
for line in genome:
	line = line.rstrip()
	if line.startswith('>'):
		if not chrom:
			chrom = line[1:]
			continue
		if chrom in psldata:
			for entry in psldata[chrom]:
				pulledseq = get_sequence(entry, seq)
				sys.stderr.write(pulledseq.upper()+'\n')
				pulledseq = revcomp(pulledseq.upper()).replace('T', 'U')
				print(entry[9], translate_seq(pulledseq), translate_seq(pulledseq[1:]), translate_seq(pulledseq[2:]))
			break
		chrom = line[1:]
		seq = ''
	else:
		seq += line
