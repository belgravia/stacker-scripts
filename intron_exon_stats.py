import sys,csv

try:
	psl = open(sys.argv[1])
	# outfilebase = sys.argv[2]
except:
	sys.stderr.write('usage: script.py psl\n')
	sys.exit()

def get_junctions(line):
	junctions = set()
	starts = [int(n) + 1 for n in line[20].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[18].split(',')[:-1]]  # for indexing pupropses
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		# junctions += [(starts[b]+sizes[b], starts[b+1])]
		junctions.add((line[13], starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_exons(line):
	blocks = set()
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]  # for indexing pupropses
	if len(starts) == 1:
		return
	for b in range(len(starts)): # block
		# junctions += [(starts[b]+sizes[b], starts[b+1])]
		blocks.add((line[13],starts[b],starts[b]+sizes[b]))
	return blocks

introns = set()
exons = set()
num_isoforms = 0
# numlines = 0
for line in psl:
	line = line.rstrip().split('\t')
	introns.update(get_junctions(line))
	exons.update(get_exons(line))
	if '-' not in line[9][-5:]:
		num_isoforms += 1
	# numlines += 1
	# if numlines % 10000 == 0:
	# 	print('hit ' + str(numlines))

sys.stderr.write('{} unique exons\n'.format(len(exons)))
sys.stderr.write('{} unique introns\n'.format(len(introns)))
sys.stderr.write('{} unique isoforms\n'.format(num_isoforms))




# with open(outfilebase + '.intron.txt' + str(level) + '.txt', 'wt') as intron_out:
# 	writer = csv.writer(intron_out, delimiter='\t')