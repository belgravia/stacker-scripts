import sys, csv

try:
	header = sys.argv[1]
except:
	print('pls specify 2 arguments, an ID (or a psl with IDs) to query in a fasta file')
	sys.exit(1)

if len(sys.argv) > 2:
	fasta = sys.argv[2]
else:
	fasta = '/scratch/alison/blat_aligned/2Dk700e_pass_fail.fa'

class FastAreader:
    def __init__(self, fname=''):
        self.fname = open(fname, 'r')
        if not fname:
            self.fname = sys.stdin
 
    def readFasta(self):
        line = self.fname.readline()
        while not line.startswith('>') :
            line = self.fname.readline()
        header = line.rstrip()
        sequence = ''
        for line in self.fname:
            if line.startswith ('>'):
                yield [header,sequence]
                header = line.rstrip()
                sequence = ''
            else :
                sequence += ''.join(line.rstrip().split()).upper()
        yield [header, sequence]


if header[-3:] == 'psl':
	psl = True
else:
	psl = False

if psl:
	print('havent actually implemented this yet')
	sys.exit(1)

fasta = FastAreader(fasta)
printed = False
for head, seq in fasta.readFasta():
	if header in head:
		print(head)
		print(seq)
		printed = True
if not printed:
	print('maybe the fasta is wrong?')
