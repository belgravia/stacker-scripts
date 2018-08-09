import sys, csv

try:
	genecountsfile = open(sys.argv[1])
	geneseqfile = open(sys.argv[2])
	base = sys.argv[3]
	elevels = 4.0  # levels of expression
	levels = [x/elevels for x in range(1, int(elevels)+1, 1)]
	size = 7  # size of kmer
except:
	sys.stderr.write('usage: script.py gene_counts gene_seq outfilenamebase \n')
	sys.exit(1)

def add_kmers(kmers, level, seq, gene, genesperlevel):
	for i in range(0, len(seq)-size):
		k = seq[i:i+size]
		if k[3].upper() != 'A': # middle fixed as an A
			continue
		if k not in kmers[level]:
			kmers[level][k] = 0
		kmers[level][k] += 1
	return kmers, genesperlevel

genecounts = []
for line in genecountsfile:
	line = line.rstrip().split('\t')
	genecounts += [(line[0], int(line[1]))]
# print(genecounts)

genecounts = sorted(genecounts, key=lambda x: x[1], reverse=True)  # sort by descending count
numgenes = float(len(genecounts))
generank = {}  # dict contains rank
rank = 0
for g in genecounts:
	generank[g[0]] = {}
	generank[g[0]]['count'] = g[1]
	generank[g[0]]['rank'] = rank
	rank += 1

kmers = dict.fromkeys(levels)
for l in kmers:
	kmers[l] = {}
	kmers[l]['kmers'] = {}
	kmers[l]['genes'] = set()
# genesperlevel = dict.fromkeys(levels)
# for l in genesperlevel:
# 	genesperlevel[l] = {}
# 	genesperlevel[l]
gene = ''
check = dict.fromkeys(levels, 0)


for line in geneseqfile:
	line = line.rstrip()
	if line.startswith('>'):
		if gene:
			level = [l for l in levels if generank[gene]['rank']/numgenes < l][0]
			check[level] += 1
			kmers = add_kmers(kmers, level, seq, gene)
		gene = line[1:]
		seq = ''
	else:
		seq += line
level = [l for l in levels if generank[gene]['rank']/numgenes < l][0]
kmers = add_kmers(kmers, level, seq)

print(check)

for level in kmers:
	with open(base + '.' + str(level) + '.txt', 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t')
		kvals = []
		for k in kmers[level]:
			writer.writerow([k, kmers[level][k]])
			kvals += [kmers[level]['kmers'][k]]
		# panel1 = sns.kdeplot(kvals, bw=0.8, color='y', shade=True)

		# panel1.set_ylabel('Density', fontsize=fs)  # 
		# # panel1.set_yticks(np.arange(0,0.5,0.1))

		# plt.savefig(base+'.'+str(level)+'.png')
for l in kmers:
	sys.stderr.write('Top {}: {} genes represented\n'.format(l, len(kmers[l]['genes']))
