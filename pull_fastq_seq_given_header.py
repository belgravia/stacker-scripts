import sys, csv

try:
	header = sys.argv[1]
	fastq = open(sys.argv[2])
except:
	sys.stderr.write('usage: script.py ID/partialID fastq\n')
	sys.exit(1)


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

n=0
found=-4
for line in fastq:
	if header in line:
		found=n
	if n < found + 4:
		print(line.rstrip())
	n += 1
		
