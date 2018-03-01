import sys, csv

try:
	psl = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py pslfile outfilename\n')
	sys.exit(1)

def overlap(coordsa, coordsb, margin=200):
	return coordsb[0] >= coordsa[0] and coordsb[0] + margin <= coordsa[1] or \
		coordsa[0] >= coordsb[0] and coordsa[0] + margin <= coordsb[1] or \
		coordsb[0] >= coordsa[0] and coordsb[1] <= coordsa[1] or \
		coordsa[0] >= coordsb[0] and coordsa[1] <= coordsb[1]

potential_chimeric = {}  # {read name: [entries]}
for line in psl:
	line = line.rstrip().split('\t')
	readname = line[9][:line[9].find('_')]
	if line[9] in potential_chimeric:
		potential_chimeric[readname] += [line]
	else:
		potential_chimeric[readname] = [line]

fusions_found = {}  # {fused genes: count}
for p in potential_chimeric:
	if len(potential_chimeric[p]) == 1:
		continue
	elif len(potential_chimeric[p]) >= 2:
		# sys.stderr.write('{} has {} alignments\n'.format(p, len(potential_chimeric[p])))
		ab = []
		for i in range(len(potential_chimeric[p])-1):
			for j in range(i+1, len(potential_chimeric[p])):
				coordi = int(potential_chimeric[p][i][11]), int(potential_chimeric[p][i][12])
				coordj = int(potential_chimeric[p][j][11]), int(potential_chimeric[p][j][12])
				if not overlap(coordi, coordj, 50):
					ab += [(i, j, coordi[1]-coordi[0]+coordj[1]-coordj[0])]
		if not ab:
			continue
		ab = sorted(ab, key=lambda x: x[2])[-1]
		i, j = ab[0], ab[1]
		coorda = int(potential_chimeric[p][i][11]), int(potential_chimeric[p][i][12])
		coordb = int(potential_chimeric[p][j][11]), int(potential_chimeric[p][j][12])
	# else:
		# coorda = int(potential_chimeric[p][0][11]), int(potential_chimeric[p][0][12])
		# coordb = int(potential_chimeric[p][1][11]), int(potential_chimeric[p][1][12])
	if overlap(coorda, coordb, 50):
		continue
	# sys.stderr.write('{} mappable regions are {} {}\n'.format(p, coorda, coordb))
	gene0 = potential_chimeric[p][0][9][potential_chimeric[p][0][9].rfind('_'):]
	gene1 = potential_chimeric[p][1][9][potential_chimeric[p][1][9].rfind('_'):]
	if gene0 == '_':
		gene0 = potential_chimeric[p][0][13]+':'+str(potential_chimeric[p][0][15])
	if gene1 == '_':
		gene1 = potential_chimeric[p][1][13]+':'+str(potential_chimeric[p][1][15])
	if coorda[1] > coordb[1]:  # coordb gene 3' end fused to 5' end of coorda gene
		fusion_name = gene0+'_'+gene1
	else:
		fusion_name = gene1+'_'+gene0
	if fusion_name not in fusions_found:
		fusions_found[fusion_name] = {}
		fusions_found[fusion_name]['count'] = 0
		fusions_found[fusion_name]['readnames'] = []
	fusions_found[fusion_name]['count'] += 1
	fusions_found[fusion_name]['readnames'] += [p]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for f in fusions_found:
		writer.writerow([f, fusions_found[f]['count'], ','.join(fusions_found[f]['readnames'])])
