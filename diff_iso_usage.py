import sys, csv
import scipy.stats as sps

try:
	psl = open(sys.argv[1])  # wt
	outfilename = sys.argv[2]
	if len(sys.argv) > 3:
		col = int(sys.argv[3])
	else:
		col = 25
	if len(sys.argv) > 4:
		matchedoutfile = sys.argv[4]
	else:
		matchedoutfile = ''
except:
	sys.stderr.write('usage: script.py psl outfilename [colname] [matchedcountoutfile]\n')
	sys.exit(1)

#ntn = {}  # name to name
#for line in matches:
#	line = line.rstrip().split('\t')
#	name1 = line[0][:line[0].find('Isoform')-1]
#	iso1 = line[0][line[0].find('Isoform'):]
#	iso1 = iso1[:iso1.find('_')]
#	# name2 = line[1][:line[1].find('Isoform')]
#	# iso2 = line[1][line[1].find('Isoform'):]
#	# iso2 = iso2[:iso2.find('_')]
#	if name1 not in ntn:
#		ntn[name1] = {}
#	ntn[name1][iso1] = line[1]

counts = {}
ever = 0

for line in psl:
	line = line.rstrip().split('\t')
	gene = line[9][line[9].rfind('_')+1:]
	if gene.count('.') == 2:
		gene = gene[:gene.rfind('.')]
	if gene not in counts:
		counts[gene] = {}
	if line[col] == 'NA':
		line[col] = 0
	if line[col+1] == 'NA':
		line[col+1] = 0
	counts[gene][line[9]] = [float(line[col]), float(line[col+1])]

if matchedoutfile:
	with open(matchedoutfile, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t')
		for gene in counts:
			for iso in counts[gene]:
				writer.writerow([gene, iso] + counts[gene][iso])

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	geneordered = sorted(counts.keys())
	for gene in geneordered:
		generes = []
		for iso in counts[gene]:
			thesecounts = counts[gene][iso]
			othercounts = [0, 0]
			for iso_ in counts[gene]:  # count up for all other isoforms of this gene
				if iso_ == iso:
					continue
				othercounts[0] += counts[gene][iso_][0]
				othercounts[1] += counts[gene][iso_][1]
			ctable = [thesecounts, othercounts]
			if ctable[0][0] + ctable[1][0] == 0 or ctable[0][1] + ctable[1][1] == 0 or not sum(ctable[1]):
				continue
			generes += [[gene, iso, sps.fisher_exact(ctable)[1]] + \
						 ctable[0] + ctable[1]]
		if not generes:
			continue
		generes = sorted(generes, key=lambda x: x[2])
		# print(generes)
		writer.writerow(generes[0])

# sys.stderr.write('overlap: {}, ever: {}\n'.format(both, ever))
