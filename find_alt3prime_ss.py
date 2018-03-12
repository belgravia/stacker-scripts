import sys, csv
import scipy.stats as sps

try:
	bed1 = open(sys.argv[1])
	bed2 = open(sys.argv[2])
	outfilename = sys.argv[3]
	fiveprimeon = True if len(sys.argv)>4 else False
except:
	print('usage: script.py junctionswt junctionsmt outfilename [5 prime]')  # all junctions (not just the ones unique to long reads)
	sys.exit()


annotmin = {}  # all annotated 3' SS on the minus strand
annotpos = {}
geneset = set()
prevgene = ''
for line in open('/pod/pstore/groups/brookslab/atang/annotation/gencode.v24.annotation.gtf'):
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand, gene = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8]
	if ty != 'exon':
		continue
	gene = gene[gene.find('gene_id')+len('gene_id')+2:]
	gene = gene[:gene.find('"')]
	geneset.add(gene)
	firstexon = not gene == prevgene
	prevgene = gene
	if strand == '+':
		if chrom not in annotpos:
			annotpos[chrom] = []
		if fiveprimeon:  # this will store the end of the last exon, unfortunately
			annotpos[chrom] += [end]  # exon ends are SS acceptor positions
			continue
		if not firstexon:
			annotpos[chrom] += [start]  # 3' end of splice site stored
	else:
		if chrom not in annotmin:
			annotmin[chrom] = []
		if not fiveprimeon:
			annotmin[chrom] += [end]
			continue
		if not firstexon:
			annotmin[chrom] += [start]


def find_wiggle(coord, annot, annot2={}, maxdist=100):
	""" Finds the distance between coordinate and the closest annotated pos in annot dict. """
	wiggle = 0
	coord = int(coord)
	while coord + wiggle not in annot and coord + wiggle not in annot2:
		if wiggle == maxdist:
			break
		if wiggle == 0:
			wiggle += 1
		elif wiggle >= 0:
			wiggle = wiggle * -1
		else:
			wiggle = (wiggle-1) * -1
	return wiggle


def bedreader(bed, junctiondict, index=0):
	for line in bed:
		line = line.rstrip().split('\t')
		if line[5] not in ['+', '-']:
			continue
		elif fiveprimeon and line[5] == '+':
			line[1], line[2] = line[2], line[1]
		elif not fiveprimeon and line[5] == '-':
			line[1], line[2] = line[2], line[1]  # reverse coordinates for junctions on - strand

		chrom, fiveprime, threeprime, name, count, strand = line
		chrom = strand + chrom
		if chrom not in junctiondict:
			junctiondict[chrom] = {}
		if fiveprime not in junctiondict[chrom]:
			junctiondict[chrom][fiveprime] = {}  # 5' end anchor
		if threeprime not in junctiondict[chrom][fiveprime]:
			junctiondict[chrom][fiveprime][threeprime] = [0,0, name]
		junctiondict[chrom][fiveprime][threeprime][index] = int(count)
	return junctiondict

alljuncs = {}
alljuncs = bedreader(bed1, alljuncs)
alljuncs = bedreader(bed2, alljuncs, 1)

ambiguous_junctions = set()
far = 0
maxdist=100
nocanonical=0
totalnumsites = 0
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in alljuncs:
		for fiveprime in alljuncs[chrom]:
			if len(alljuncs[chrom][fiveprime]) == 1:  # if there is only one 3' end, not an alt 3' junction
				continue
			entries = []
			for threeprime1 in alljuncs[chrom][fiveprime]:  # for each potential cryptic ss
				if sum(alljuncs[chrom][fiveprime][threeprime1][:2]) < 10 :  # if this ss is not really used
					continue
				allothercounts = [0,0]
				oro3p = ('', 0)  # overrepresented other 3 prime site, used to calculate distance from cryptic SS
				ambi = False
				farused = False
				for threeprime2 in alljuncs[chrom][fiveprime]:
					if threeprime1 == threeprime2:
						continue
					if abs(int(threeprime1) - int(threeprime2)) < 9:  # if the ss are close together -> ambiguity
						ambiguous_junctions.add(chrom+':'+threeprime1)
						ambi = True
						break
					if abs(int(threeprime1) - int(threeprime2)) > 200:  # likely an exon skipping and not a alt 3/5' site
						if not farused:
							far += 1
							farused = True
						continue
					if int(alljuncs[chrom][fiveprime][threeprime2][0]) > oro3p[1]:
						oro3p = (threeprime2, int(alljuncs[chrom][fiveprime][threeprime2][0]))
					allothercounts[0] += alljuncs[chrom][fiveprime][threeprime2][0]
					allothercounts[1] += alljuncs[chrom][fiveprime][threeprime2][1]
				if sum(allothercounts) <= 5 or ambi:  # if not enough alternative ss usage
					continue
				ctable = [alljuncs[chrom][fiveprime][threeprime1][:2], allothercounts]
				name = alljuncs[chrom][fiveprime][threeprime1][2]
				if oro3p[1] == 0:
					continue
				# if ctable[0][1] > ctable[0][0]:  # this junction has more counts in the mutant
				totalnumsites += 1
				if chrom[0] == '-':
					closest_annotated = find_wiggle(threeprime1, annotmin, maxdist=150)
					if closest_annotated == 150:
						nocanonical+=1
						diff = min(closest_annotated, int(oro3p[0])-int(threeprime1))

					entries += [[chrom[1:], threeprime1, fiveprime, sps.fisher_exact(ctable)[1], 
					chrom[0]] + ctable[0] + ctable[1] + [name] + [diff] + [oro3p[0]]]
					# writer.writerow([chrom[1:], threprime1, fiveprime, sps.fisher_exact(ctable)[1], 
					# chrom[0]] + ctable[0] + ctable[1] + [name] + [int(oro3p[0])-int(threeprime1)] + [oro3p[0]])
				else:
					closest_annotated = find_wiggle(threeprime1, annotpos, maxdist=150)
					if closest_annotated == 150:
						nocanonical+=1
						diff = min(closest_annotated, int(threeprime1)-int(oro3p[0]))	
					entries += [[chrom[1:], fiveprime, threeprime1, sps.fisher_exact(ctable)[1], 
					chrom[0]] + ctable[0] + ctable[1] + [name] + [diff] + [oro3p[0]]]
					# writer.writerow([chrom[1:], fiveprime, threeprime1, sps.fisher_exact(ctable)[1], 
					# chrom[0]] + ctable[0] + ctable[1] + [name] + [int(threeprime1)-int(oro3p[0])] + [oro3p[0]])
			if entries:
				entries = sorted(entries, key=lambda x: x[3])
				writer.writerow(entries[0])
print(len(ambiguous_junctions))
print(far)
print(nocanonical, totalnumsites)