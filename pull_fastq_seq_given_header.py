import sys, csv

try:
	header = sys.argv[1]
	fastq = open(sys.argv[2])
except:
	sys.stderr.write('usage: script.py ID/partialID/psl fastq > output\n')
	sys.exit(1)


headers_keep = []
if header[-3:] == 'psl':
    psl = True
    for line in open(header):
        line = line.rstrip().split('\t')
        headers_keep += [line[9]]
else:
    psl = False

n = 0
if not psl:
    found=-4
    for line in fastq:
    	if header in line:
    		found=n
    	if n < found + 4:
    		print(line.rstrip())
    	n += 1
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

