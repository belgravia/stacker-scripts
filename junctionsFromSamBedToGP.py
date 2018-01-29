import sys, csv

try:
	bed = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py bed outgp\n')
	sys.exit(1)

exondist = {}  # distribution of the number of exons
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in bed:
		line = line.rstrip().split('\t')
		relblockstarts = [int(x) for x in line[11].split(',')]
		blocksizes = [int(x) for x in line[10].split(',')]
		blockstarts = [x + int(line[1]) for x in relblockstarts]
		blockends = ','.join([str(x+y) for (x,y) in zip(blockstarts, blocksizes)])+','
		blockstarts = ','.join([str(x) for x in blockstarts])+','
		# if int(line[4]) not in exondist:
		# 	exondist[int(line[4])] = 1
		# else:
		# 	exondist[int(line[4])] += 1
		writer.writerow([line[3], line[0], line[5], line[1], line[2], line[1], line[2], line[9], 
			blockstarts, blockends, line[4]])

# print(exondist)