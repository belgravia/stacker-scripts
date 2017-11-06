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
		blockstarts = [int(x) for x in line[20].split(',')[:-1]]
		blocksizes = [int(x) for x in line[18].split(',')[:-1]]
		newstarts = []
		newsizes = []
		for i in range(len(blocksizes)):
			size = blocksizes[i]
			if size >= minsize:
				newsizes += [size]
				newstarts += [blockstarts[i]]
		line[20] = ','.join([str(s) for s in newstarts]) + ','
		line[18] = ','.join([str(s) for s in newsizes]) + ','
		writer.writerow(line)

