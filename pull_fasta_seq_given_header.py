import sys, csv

try:
    header = sys.argv[1]
    fasta = sys.argv[2]
except:
    sys.stderr.write('script.py header/psl fasta\n')
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
            else :
                sequence += ''.join(line.rstrip().split()).upper()
        yield [header, sequence]


headers_keep = []
if header[-3:] == 'psl':
    psl = True
    for line in open(header):
        line = line.rstrip().split('\t')
        headers_keep += [line[9]]
else:
    psl = False

fasta = FastAreader(fasta)
for head, seq in fasta.readFasta():
    if not psl and header in head:
        print(head)
        print(seq)
    elif psl and head[1:] in headers_keep:
        print(head)
        print(seq)
    # sys.exit()



