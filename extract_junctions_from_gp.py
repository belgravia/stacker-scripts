import sys, csv

try:
	gp = open(sys.argv[1]) 
	outfilename = sys.argv[2]
except:
	print('usage: script.py corrected.gp outfilename')
	print('extracts junctions from a gp file, written to complement jeltje script')
	sys.exit(1)

counts = {}

counts_unique = {}

doesthishappen = 0
for line in gp:
	line = line.rstrip().split('\t')
	blockstarts = line[8].split(',')[:-1]
	blockends = line[9].split(',')[:-1]
	for s in range(len(blockstarts)-1):
		jend = blockstarts[s+1] #+ blocksizes[s]
		jstart = blockends[s] #blockstarts[s+1]
		junction = (line[1], int(jstart), int(jend), line[2])  # chr #jstart #jend #strand
		# junction = (line[1], int(jstart), int(jend), '+')  # chr #jstart #jend #strand
		if junction not in counts:
			counts[junction] = 1
		else:
			counts[junction] += 1

		if junction not in counts_unique:
			counts_unique[junction] = [1, set()]  # read name
			counts_unique[junction][1].add(line[0])
		else:
			if line[0] in counts_unique[junction][1]:
				doesthishappen += 1
				continue 
			else:
				counts_unique[junction][0] += 1
				counts_unique[junction][1].add(line[0])

junctions = counts_unique.keys()
junctions = sorted(junctions)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for j in junctions:
		writer.writerow([j[0], str(j[1]), str(j[2]), j[0]+':'+str(j[1])+'-'+str(j[2]), counts_unique[j][0], j[3]])




