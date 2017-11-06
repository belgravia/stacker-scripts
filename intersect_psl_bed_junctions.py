import sys, csv

try:
	psl = open(sys.argv[1])
	bed = open(sys.argv[2])  # chr pos pos 
	outfilename = sys.argv[3]
except:
	print('usage: script.py pslfile bedfile outfilename [acceptable_splice_window]')
	print('quantifies the number of reads in the psl that support splice junctions specified in the bedfile')
	sys.exit(1)

try:
	err = int(sys.argv[4]) + 1
except:
	err = 6

print(err)
diffs = dict.fromkeys(range(-10, 10), 0)
def spans_junction(starts, blocksizes, junction_0, junction_1):  # junction_N are intronic positions defining a junction
	first_exon_coverage = second_exon_coverage = intron_coverage = False
	if starts[0] < junction_0:
		for i in range(len(starts)):
			d = junction_0 - (starts[i] + blocksizes[i])
			if abs(d) < 10:
				diffs[d] += 1
			first_exon_coverage = first_exon_coverage or abs(junction_0 - (starts[i] + blocksizes[i])) in range(err)
			second_exon_coverage = second_exon_coverage or abs(junction_1 - starts[i]) in range(err)
			intron_coverage = intron_coverage or (starts[i] < junction_1 and starts[i] > junction_0)
	return first_exon_coverage and second_exon_coverage and not intron_coverage

junctions = {}
junction_support = {}
for line in bed:
	line = line.rstrip().split('\t')
	if line[0] not in junctions:  # chr
		junctions[line[0]] = []
	junctions[line[0]] += [(line[3], int(line[1]), int(line[2]))]  # chrom, genename, coord0, coord1
	if line[3] not in junction_support:  # gene
		junction_support[line[3]] = {}
	junction_support[line[3]][(int(line[1]), int(line[2]))] = 0  # initialize the counter for each gene

for chrom in junctions:
	junctions[chrom].sort(key=lambda x: x[1])

e = 1
for line in psl:
	line = line.rstrip().split('\t')
	if not e % 5e4:
		print(str(e) + ' PSL entries processed')
	e += 1
	chrom = line[13]
	if chrom not in junctions:  # if it is not a relevant chromosome
		continue
	read_start = int(line[15])
	read_end = int(line[16])
	tstarts = line[20].split(',')[:-1]
	blocksizes = line[18].split(',')[:-1]
	tstarts = [int(s) for s in tstarts]
	blocksizes = [int(s) for s in blocksizes]
	# for i in range(len(blocksizes)):
		# read_alignments += [(tstarts[i], tstarts[i]+blocksizes[i])]
	# read_alignments += [(read_start, read_end)]
	# junc = junctions[chrom][len(junctions[chrom])/2]
	# if junc[1] < read_start:  # give some uncertainty here
		# if spans_junction(junctions[chrom][len(junctions[chrom])/2-10:]):  # use second half 
	# elif junc[1] == read_start:  # use 

	# elif junc[2] > read_end:  # use first half
		# if 

	for j in junctions[chrom]:   # insert chrom here
		if spans_junction(tstarts, blocksizes, j[1], j[2]):
			junction_support[j[0]][(j[1], j[2])] += 1

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for gene in junction_support:
		for junction in junction_support[gene]:
			writer.writerow([gene, junction[0], junction[1], junction_support[gene][junction]])
print(str(e) + ' PSL entries')
print(diffs)