import sys, csv

try:
	psl = open(sys.argv[1])
	maxnum = int(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py psl tes/tss_threshold outfilename \n')
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
if '_0.' in sys.argv[1]:
	psl.readline()

n = 0
for line in psl:
	line = tuple(line.rstrip().split('\t'))
	chrom = line[13]
	if chrom not in isoforms:
		isoforms[chrom] = {}
	junctions = get_junctions(line)
	tss, tes = get_start_end(line)
	counts = 0 if len(line) < 36 else int(line[35])
	for i in isoforms[chrom]:  # does a similar isoform already exist
		if junctions == isoforms[chrom][i]['junctions']:  # must have same junctions
			if abs(tss - isoforms[chrom][i]['start']) < maxnum and \
				abs(tes - isoforms[chrom][i]['end']) < maxnum:  # tes and tss within maxnum of another isoform are considered the same
				tss = min(tss, isoforms[chrom][i]['start'])
				tes = max(tes, isoforms[chrom][i]['end'])
				counts = isoforms[chrom][i]['counts']
				isoforms[chrom].pop(i)
				break
	isoforms[chrom][line] = {}
	isoforms[chrom][line]['junctions'] = junctions
	isoforms[chrom][line]['start'] = tss
	isoforms[chrom][line]['end'] = tes
	isoforms[chrom][line]['counts'] = 1 + counts
	n += 1
	if n % 10000 == 0:
		print(n)
# print(isoforms)
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in isoforms:
		for iso in isoforms[chrom]:
			writer.writerow(list(iso) + [isoforms[chrom][iso]['counts']])
			# if isoforms[chrom][iso]['counts'] > minnum:
				# writer.writerow(iso)
