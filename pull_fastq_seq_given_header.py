import sys, csv

try:
    header = sys.argv[1]
    sam = sys.argv[2][-4:] == '.sam'
    fastq = open(sys.argv[2])
except:
	sys.stderr.write('usage: script.py ID/partialID/psl fastq/sam > output\n')
	sys.exit(1)


headers_keep = []
psl = header[-4:] == '.psl'

n = 0
if not psl:
    found =- 4
    sys.stderr.write('searching for {} in headers of fastq file\n'.format(header))
    for line in fastq:
        if header in line:
            found = n
        if n < found + 4:
    		print(line.rstrip())
    	n += 1
    sys.exit()

if sam:
    readinfo = {}
    seen = set()
    for line in fastq:  # samfile
        line = line.rstrip().split()
        qname, flag, tname, pos, cigar, seq, qual = line[0], int(line[1]), line[2], int(line[3]), line[5], line[9], line[10]
        readinfo[qname] = seq

    for line in open(header):
        line = line.rstrip().split('\t')
        name = line[9]
        if name in seen:
            continue
        seen.add(name)
        print('>' + name)
        print(readinfo[name])

for line in open(header):  # psl file
    line = line.rstrip().split('\t')
    headers_keep += [line[9]]



    sys.exit()


found = False
for line in fastq:
    if n % 4 == 0:
        line = line.rstrip().split()[0]
        if line.startswith('@'):
            line = line[1:]
        if line in headers_keep:
            found = True
            # print(line.rstrip())
            print('>' + line.rstrip())
    elif n % 4 == 1 and found:
        print(line.rstrip())   
    elif n % 4 == 2 and found:
        # print(line.rstrip())  # uncomment if fastq output is desired
        pass
    elif found:
        # print(line.rstrip())
        found = False
    n += 1

