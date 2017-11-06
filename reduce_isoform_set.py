import sys, csv

try:
	psl = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py pslfile outfilename\n')
	sys.exit(1)

def get_junctions(line):
	junctions = set()
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	for b in range(len(starts)-1): # block
		# junctions += [(starts[b]+sizes[b], starts[b+1])]
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_start_end(line):
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	return starts[0], starts[-1]+sizes[-1]

isoforms = {}
for line in psl:
	line = tuple(line.rstrip().split('\t'))
	chrom = line[13]
	if chrom not in isoforms:
		isoforms[chrom] = {}
	junctions = get_junctions(line)
	tss, tes = get_start_end(line)
	for i in isoforms[chrom]:
		if junctions == isoforms[chrom][i]['junctions']:
			if abs(tss - isoforms[chrom][i]['start']) < 15 and \
				abs(tes - isoforms[chrom][i]['end']) < 15:
				tss = min(tss, isoforms[chrom][i]['start'])
				tes = max(tes, isoforms[chrom][i]['end'])
				isoforms[chrom].pop(i)
				break
	isoforms[chrom][line] = {}
	isoforms[chrom][line]['junctions'] = junctions
	isoforms[chrom][line]['start'] = tss
	isoforms[chrom][line]['end'] = tes

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in isoforms:
		for i in isoforms[chrom]:
			writer.writerow(list(i))