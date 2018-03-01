import sys, csv

try:
    psl = open(sys.argv[1])
    wiggle = int(sys.argv[2])
    outfilename = sys.argv[3]
except:
    sys.stderr.write('usage: script.py pslfile outfilename\n')
    sys.exit(1)

with open(outfilename, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    for line in psl:
        line = line.rstrip().split('\t')
        tlen, alignend = int(line[14]), int(line[16])
        if tlen - alignend < wiggle:
        	writer.writerow(line)
