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
	if len(starts) == 1:
		return
	for b in range(len(starts)-1): # block
		# junctions += [(starts[b]+sizes[b], starts[b+1])]
		junctions.add((starts[b]+sizes[b], starts[b+1]))
	return junctions

def get_start_end(line):
	starts = [int(n) for n in line[20].split(',')[:-1]]
	sizes = [int(n) for n in line[18].split(',')[:-1]]
	return starts[0], starts[-1]+sizes[-1]

def find_best_site(sites):
	total = float(sum(list(sites.values())))
	if total < 3:
		return ''
	nearby = dict.fromkeys(sites, 0)  # site, distance away
	bestsite = (0, 0)
	for s in sites:  # calculate number of reads supporting this site within maxnum
		for s_ in sites:
			if abs(s - s_) <= maxnum:
				nearby[s] += sites[s_]
		if nearby[s] > bestsite[1]:
			bestsite = (s, nearby[s])
	return [bestsite]  # no alternative sites allowed. sorry :c 
	if total - bestsite[1] < 3:  # best
		return [bestsite]
	for s in sites:  # remove bestsite reads
		if abs(s - bestsite[0]) <= maxnum:
			nearby[s] = 0
	other_sites = [bestsite]
	for s in sites:  # find other sites with sufficient coverage
		close = False
		if nearby[s] < 3:
			continue
		worseothersite = ''
		for os in other_sites:
			if abs(os[0] - s) <= maxnum * 1.5:
				if nearby[s] > os[1]:
					worseothersite = os
				close = True
				break
		if worseothersite:
			other_sites.remove(worseothersite)
			other_sites += [(s, nearby[s])]  # replacing 
		if not close:
			other_sites += [(s, nearby[s])]  # adding
	other_sites = sorted(other_sites, key=lambda x: x[1])[-2:]
	return other_sites

def edit_line(line, tss, tes, blocksize=''):
	if blocksize:
		line[18] = str(blocksize) + ','
		line[20] = str(tss) + ','
		line[15] = tss
		line[16] = tes
		return line

	bsizes = [int(x) for x in line[18].split(',')[:-1]]
	bstarts = [int(x) for x in line[20].split(',')[:-1]]
	tstart = int(line[15])  # current chrom start
	tend = int(line[16])
	bsizes[0] += tstart - tss
	bsizes[-1] += tes - tend
	if bsizes[0] < 0:
		print(tss, tstart, tend, line)
		return ''
	bstarts[0] = tss
	line[15] = tss
	line[16] = tes
	line[18] = ','.join([str(x) for x in bsizes])+','
	line[20] = ','.join([str(x) for x in bstarts])+','
	return line

def single_exon_pairs(sedict):  # incomplete
	tss_nearby = dict.fromkeys(sedict['tss'], 0)
	tes_nearby = dict.fromkeys(sedict['tes'], 0)
	for tss in sedict['tss']:
		for tss_ in sedict['tss']:
			tss_nearby[tss] += sedict['tss'][tss_]
	for tes in sedict['tes']:
		for tes_ in sedict['tes']:
			tes_nearby[tes] += sedict['tes'][tes_]
	alltss = sorted(tss_nearby.keys())
	alltes = sorted(tes_nearby.keys())
	tssi = 0
	tesi = 0
	tss_stack = []
	pairs = []
	while tssi < len(alltss) and tesi < len(alltes):
		if tes_nearby[alltes[tesi]] < 3:
			tssi += 1
			continue
		if alltes[tesi] < alltss[tssi]:
			pairs += []
			tss_stack = [alltss[tssi]]
		else:
			tss_stack += [alltss[tssi]]
	return pairs

isoforms = {}
n = 0
singleexon = {}
for line in psl:
	line = tuple(line.rstrip().split('\t'))
	chrom = line[13]
	tss, tes = get_start_end(line)
	junctions = get_junctions(line)

	if not junctions:
		junctions = 'se'  # all single exons will be considered together
		if chrom not in singleexon:
			singleexon[chrom] = {}
			singleexon[chrom]['tss'] = {}
			singleexon[chrom]['tes'] = {}
		if tss not in singleexon[chrom]['tss']:
			singleexon[chrom]['tss'][tss] = 0
		singleexon[chrom]['tss'][tss] += 1
		if tes not in singleexon[chrom]['tes']:
			singleexon[chrom]['tes'][tes] = 0
		singleexon[chrom]['tes'][tes] += 1
		continue

	junctions = str(sorted(list(junctions)))  # hashable but still unique
	if chrom not in isoforms:
		isoforms[chrom] = {}
	if junctions not in isoforms[chrom]:
		isoforms[chrom][junctions] = {}
		isoforms[chrom][junctions]['tss'] = {}
		isoforms[chrom][junctions]['tes'] = {}
		isoforms[chrom][junctions]['line'] = line
	if tss not in isoforms[chrom][junctions]['tss']:
		isoforms[chrom][junctions]['tss'][tss] = 0
	isoforms[chrom][junctions]['tss'][tss] += 1
	if tes not in isoforms[chrom][junctions]['tes']:
		isoforms[chrom][junctions]['tes'][tes] = 0
	isoforms[chrom][junctions]['tes'][tes] += 1

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in isoforms:
		for jset in isoforms[chrom]:
			tss = find_best_site(isoforms[chrom][jset]['tss'])
			tes = find_best_site(isoforms[chrom][jset]['tes'])
			if not tes:
				continue
			line = isoforms[chrom][jset]['line']
			# if len(tss[0]) and len(tes[0]) == 1:
			# 	line = edit_line(line, tss[0], tes[0])
			# 	writer.writerow(line + [max(tss[1], tes[1])])
			# 	if tss[1] != tes[1]:
					# print(tss[1], tes[1],isoforms[chrom][jset]['tss'], isoforms[chrom][jset]['tes'] )
			# else:
			if len(tes) > 1 and len(tss) > 1:
				print(line, tes, tss)
			for tes_ in tes:
				for tss_ in tss:
					if tes_[0] < 0 or tss_[0] < 0:
						print(tes_, tss_)
					templine = edit_line(list(line), tss_[0], tes_[0])
					if not templine:
						continue
					writer.writerow(templine + [max(tss_[1], tes_[1])])
	# for chrom in singleexon:
	# 	pairs = single_exon_pairs(singleexon[chrom])
	# 	for p in pairs:
	# 		tss, tes, freq = p
	# 		templine = edit_line(list(line), tss, tes, tes-tss)
	# 		if not templine:
	# 			continue
	# 		writer.writerow(templine + [freq])







