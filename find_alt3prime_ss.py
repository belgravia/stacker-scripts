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

def bedreader(bed, junctiondict, index=0):
	for line in bed:
		line = line.rstrip().split('\t')
		if line[5] not in ['+', '-']:
			continue
		elif fiveprimeon:
			if line[5] == '+':
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

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in alljuncs:
		for fiveprime in alljuncs[chrom]:
			if len(alljuncs[chrom][fiveprime]) == 1:  # if there is only one 3' end, not an alt 3' junction
				continue
			for threeprime1 in alljuncs[chrom][fiveprime]:
				if sum(alljuncs[chrom][fiveprime][threeprime1][:2]) <= 10:
					continue
				allothercounts = [0,0]
				oro3p = ('', 0)  # overrepresented other 3 prime site, used to calculate distance from cryptic SS and canonical
				for threeprime2 in alljuncs[chrom][fiveprime]:
					if threeprime1 == threeprime2:
						continue
					if int(alljuncs[chrom][fiveprime][threeprime2][0]) > oro3p[1]:
						oro3p = (threeprime2, int(alljuncs[chrom][fiveprime][threeprime2][0]))
					allothercounts[0] += alljuncs[chrom][fiveprime][threeprime2][0]
					allothercounts[1] += alljuncs[chrom][fiveprime][threeprime2][1]
				if sum(allothercounts) <= 1:
					continue
				ctable = [alljuncs[chrom][fiveprime][threeprime1][:2], allothercounts]
				name = alljuncs[chrom][fiveprime][threeprime1][2]
				if oro3p[1] == 0:
					continue
				# if ctable[0][1] > ctable[0][0]:  # this junction has more counts in the mutant
				if chrom[0] == '-':
					writer.writerow([chrom[1:], fiveprime, threeprime1, sps.fisher_exact(ctable)[1], 
					chrom[0]] + ctable[0] + ctable[1] + [name] + [int(oro3p[0])-int(threeprime1)] + [oro3p[0]])
				else:
					writer.writerow([chrom[1:], fiveprime, threeprime1, sps.fisher_exact(ctable)[1], 
					chrom[0]] + ctable[0] + ctable[1] + [name] + [int(threeprime1)-int(oro3p[0])] + [oro3p[0]])

