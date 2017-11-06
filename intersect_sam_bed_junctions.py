import sys, csv

try:
	sam = open(sys.argv[1])
	bed = open(sys.argv[2])  # chr pos pos 
	outfilename = sys.argv[3]
except:
	print('usage: script.py samfile bedfile outfilename')
	sys.exit(1)

bag = set()
closerbag = set()
diffs = dict.fromkeys(range(-10, 10), 0)
def spans_junction(read_start, cigar, junction_0, junction_1):  # junction_N are intronic positions defining a junction
	if read_start > junction_0:
		return False
	if cigar.count('M') != 2: #cigar.count('M') > 2 or cigar.count('M') <= 1: 
		bag.add(cigar)
		return False
	# add a try statement here
	try:
		matchlen_1 = int(cigar[:cigar.find('M')])
		intronlen = int(cigar[cigar.find('M')+1:cigar.find('N')])
		matchlen_2 = int(cigar[cigar.find('N')+1:cigar.rfind('M')])  # not really necessary
		read_0 = read_start + matchlen_1  # exon junction coordinates
		read_1 = read_0 + intronlen
	except ValueError:
		closerbag.add(cigar)
		return False
	d = junction_0 - (read_0)
	if abs(d) < 10:
		diffs[d] += 1	
	return abs(junction_0 - (read_0)) in range(2) and abs(junction_1 - (read_1)) in range(2)

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

i = 0
for line in sam:
	if line.startswith('@'):
		continue
	line = line.rstrip().split('\t')
	if not i % 1e6:
		print(str(i) + ' SAM entries processed')
	i += 1
	if line[2] not in junctions:
		continue
	read_start = int(line[3])
	cigar = line[5]
	upper = len(junctions[line[2]])
	pos = upper//2
	lower = 0
	while True:
		midjunc = junctions[line[2]][pos]
		# print(midjunc[1], read_start, pos, lower, upper)
		if abs(upper - lower) < 6:
			break
		if midjunc[1]-1 < read_start:
			lower = pos
			pos = (upper + pos) // 2
		elif midjunc[1] + 1 > read_start + 100:
			upper = pos
			pos = (lower + pos) // 2
		else:
			break
	relevant_juncs = junctions[line[2]][pos-3:pos+3]
	# if midjunc[1] - 1 < read_start:  # cut in half search
	# 	relevant_juncs = junctions[line[2]][len(junctions[line[2]])//2-10:]
	# else:
	# 	relevant_juncs = junctions[line[2]][:len(junctions[line[2]])//2+10]
	for j in relevant_juncs:   # insert chrom here
		if spans_junction(read_start, cigar, j[1], j[2]):
			junction_support[j[0]][(j[1], j[2])] += 1

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for gene in junction_support:
		for junction in junction_support[gene]:
			writer.writerow([gene, junction[0], junction[1], junction_support[gene][junction]])
print(bag)
print('--------')
print(closerbag)
print(diffs)
