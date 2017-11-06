import sys, csv

try:
	bed = open(sys.argv[1])
	outfilename1 = sys.argv[2]
	outfilename2 = sys.argv[3]
except:
	print('usage: script.py corrected_bedfile outfilename5 outfilename3')
	sys.exit(1)

fiveprime = {}
fiveprimerev = {}
for line in bed:
	line = line.rstrip().split('\t')
	# find junctions and add to fiveprime dictionary
	strand = line[5]
	chrom = line[0]
	if strand == '+':
		five = line[1]
		three = line[2]
		if chrom not in fiveprime:
			fiveprime[chrom] = {}
		if five not in fiveprime[chrom]:
			fiveprime[chrom][five] = {}
		if three not in fiveprime[chrom][five]:
			fiveprime[chrom][five][three] = 0
		fiveprime[chrom][five][three] = line[4]
	elif strand == '-':
		five = line[2]
		three = line[1]
		if chrom not in fiveprimerev:
			fiveprimerev[chrom] = {}
		if five not in fiveprimerev[chrom]:
			fiveprimerev[chrom][five] = {}
		if three not in fiveprimerev[chrom][five]:
			fiveprimerev[chrom][five][three] = 0
		fiveprimerev[chrom][five][three] = line[4]

dist = {}
with open(outfilename1, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in fiveprime:
		for five in fiveprime[chrom]:
			if len(fiveprime[chrom][five]) > 1:
				for three in fiveprime[chrom][five]:
					writer.writerow([chrom, five, three, 
						fiveprime[chrom][five][three]])
					num = int(fiveprime[chrom][five][three])
					if num not in dist:
						dist[num] = 1
					else:
						dist[num] += 1
with open(outfilename2, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in fiveprimerev:
		for five in fiveprimerev[chrom]:
			if len(fiveprimerev[chrom][five]) > 1:
				for three in fiveprimerev[chrom][five]:
					writer.writerow([chrom, three, five, 
						fiveprimerev[chrom][five][three]])
					a = int(fiveprimerev[chrom][five][three])
					num = int(fiveprimerev[chrom][five][three])
					if num not in dist:
						dist[num] = 1
					else:
						dist[num] += 1

print(dist)  # curious about how often junctions are being counted

