import sys
class FastAreader:
    ''' This class taken from BME 160 Lab 4 starter code.
    Provides reading of a FastA file.
    object attribute fname: standard input
    methods readFasta(): returns header and sequence as strings.
    Author: David Bernick, Modifications: me
    Date: April 19, 2013 '''
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
        if not fname:
            self.fname = sys.stdin
 
    def readFasta(self):
        ''' returns each FastA record as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header. '''
        line = self.fname.readline()
        while not line.startswith('>') :
            line = self.fname.readline()
        header = line[1:].rstrip()
        sequence = ''
        for line in self.fname:
            if line.startswith ('>'):
                yield [header,sequence]
                header = line[1:].rstrip()
                sequence = ''
            else :
                sequence += ''.join(line.rstrip().split()).upper()
        yield [header, sequence]

if len(sys.argv) > 1:
	fname = sys.argv[1]
	reader = FastAreader(fname)
else:
	reader = FastAreader()

sizes = []
for header, sequence in reader.readFasta():
	sizes += [len(sequence)]

for s in sizes:
	print(s) 







