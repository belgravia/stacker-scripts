import sys,csv

try:
	psl = open(sys.argv[1])  # psl file of all isoforms from FLAIR
except:
	sys.stderr.write('usage: script.py psl\n')
	sys.exit()

def get_junctions(line):
	junctions = set()
	starts = [int(n) + 1 for n in line[20].split(',')[:-1]]
	sizes = [int(n) - 1 for n in line[18].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		junctions.add((line[13], starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_exons(line):
	blocks = set()
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	if len(starts) == 1:
		return
	for b in range(len(starts)):
		blocks.add((line[13],starts[b],starts[b]+sizes[b]))
	return blocks

introns = set()  # set of all unique introns found in isoform file
exons = set()  # set of all unique exons found in isoform file
num_isoforms = 0  # number of isoforms with a unique set of splice junctions used
for line in psl:
	line = line.rstrip().split('\t')
	introns.update(get_junctions(line))
	exons.update(get_exons(line))
	# if '-' not in line[9][-5:]:  # - indicates redundant redundant splice junctions
	num_isoforms += 1  # consider same splice junction isiforms with different TSS/TES as different 

sys.stderr.write('{} unique exons\n'.format(len(exons)))
sys.stderr.write('{} unique introns\n'.format(len(introns)))
sys.stderr.write('{} unique isoforms\n'.format(num_isoforms))