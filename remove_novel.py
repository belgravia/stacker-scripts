import sys

try:
    gp = open(sys.argv[1])
    psl = open(sys.argv[2])
except:
    sys.stderr.write('script.py junctions.gp psl > nonovels.psl \n')
    sys.exit(1)

annotated = {}
for line in gp:
	line = line.rstrip().split('\t')
	# chrom, coord = line[0].split(':')
	chrom = line[1]
	if chrom not in annotated:
		annotated[chrom] = set()
	blockstarts = [int(n) for n in line[8].split(',')[:-1]][1:]
	blockends = [int(n) for n in line[9].split(',')[:-1]][:-1]

	# print(line[8].split(',')[-1], blockstarts, blockends)
	# sys.exit()
	for start, end in zip(blockstarts, blockends):
		annotated[chrom].add((end, start))
	# coord1, coord2 = [int(n) for n in coord.split('-')]
	# annotated[chrom].add((coord1, coord2))


# print(annotated)
lastchrom = ''
lastjunc = ''
notfound = set()
for line in psl:
	line = line.rstrip().split('\t')
	# if (str(int(line[1])+1), line[2]) not in annotated[line[0]]:
	# 	continue
	chrom = line[13]
	strand = line[8]
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]		
	validjuncs = True
	for b in range(len(starts)-1):
		junction = (starts[b]+sizes[b], starts[b+1])  # lots of indexing errors
		junction2 = (junction[0]-1, junction[1])
		junction3 = (junction[0]+1, junction[1])
		# junction3 = ()
		# print(junction, chrom)
		# sys.exit()
		if junction not in annotated[chrom] and junction2 not in annotated[chrom] and junction3 not in annotated[chrom]:
			validjuncs = False
			lastchrom = chrom
			lastjunc = junction3
			notfound.add((junction, chrom, strand))
			break
	if validjuncs:
		print('\t'.join(line))
		# pass
# 108201923, 108202982
# if notfound:
# 	sys.stderr.write('not found:\n')
# for n in notfound:
# 	sys.stderr.write('{}\n'.format(n))





