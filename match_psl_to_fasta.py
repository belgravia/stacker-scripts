#hello
import sys, csv

class FastAreader:
    ''' This class taken from BME 160 Lab 4 starter code.
    Provides reading of a FastA file.
    object attribute fname: standard input
    methods readFasta(): returns header and sequence as strings.
    Author: David Bernick, Modifications: me
    Date: April 19, 2013 '''
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = open(fname, 'r')
        if not fname:
            self.fname = sys.stdin
 
    def readFasta(self):
        ''' returns each FastA record as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header. '''
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

try:
    psl = sys.argv[1]
    fasta = sys.argv[2]
except:
    print('needs a psl and fasta, will output a fasta with only the reads in the psl')
    sys.exit(1)

pslnames = set()
for line in open(psl, 'r'):
    line = line.rstrip().split('\t')
    pslnames.add(line[9])

# namesseen = set()
fasta = FastAreader(fasta)
with open(sys.argv[1][:-4] + '.fa', 'wt') as fastafile:
    writer = csv.writer(fastafile, delimiter='\t')
    for head, seq in fasta.readFasta():
        head = head.split()[0]

        if head[1:] in pslnames:
            writer.writerow([head[:head.find('_')]])
            writer.writerow([seq])
            pslnames = pslnames - set([head[1:]])

if len(pslnames):
    print('names did not see: ')
    if len(pslnames) > 10:
        print(pslnames)
    else:
        print(pslnames)
