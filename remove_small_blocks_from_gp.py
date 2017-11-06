import sys, csv

try:
	psl = open(sys.argv[1])
	outfilename = sys.argv[2]
	if len(sys.argv) > 3:
		minsize = int(sys.argv[3])
	else:
		minsize = 30
except:
	sys.stderr.write('usage: script.py pslfile outfilename [minsize=30] \n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		blockstarts = [int(x) for x in line[8].split(',')[:-1]]
		blockends = [int(x) for x in line[9].split(',')[:-1]]
		newstarts, newends = [], []
		for i in range(len(blockstarts)):
			size = blockends[i] - blockstarts[i]
			if size >= minsize:
				newstarts += [blockstarts[i]]
				newends += [blockends[i]]
		line[8] = ','.join([str(s) for s in newstarts]) + ','
		line[9] = ','.join([str(s) for s in newends]) + ','
		writer.writerow(line)

