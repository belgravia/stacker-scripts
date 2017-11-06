import sys, csv

try:
	bed = open(sys.argv[1])
	outfilename1 = sys.argv[2]
	outfilename2 = sys.argv[3]
except:
	print('usage: script.py corrected_bedfile outfilename5 outfilename3')
	print('same code as find_3prime_ss.py except "anchors" on 3 prime ends ')
	sys.exit(1)

threeprime = {}
threeprimerev = {}
for line in bed:
	line = line.rstrip().split('\t')
	# find junctions and add to threeprime dictionary
	strand = line[5]
	chrom = line[0]
	if strand == '+':
		five = line[1]
		three = line[2]
		if chrom not in threeprime:
			threeprime[chrom] = {}
		if three not in threeprime[chrom]:
			threeprime[chrom][three] = {}
		if five not in threeprime[chrom][three]:
			threeprime[chrom][three][five] = 0
		threeprime[chrom][three][five] = line[4]
	elif strand == '-':
		five = line[2]
		three = line[1]
		if chrom not in threeprimerev:
			threeprimerev[chrom] = {}
		if three not in threeprimerev[chrom]:
			threeprimerev[chrom][three] = {}
		if five not in threeprimerev[chrom][three]:
			threeprimerev[chrom][three][five] = 0
		threeprimerev[chrom][three][five] = line[4]

dist = {}
with open(outfilename1, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in threeprime:
		for three in threeprime[chrom]:
			if len(threeprime[chrom][three]) > 1:
				for five in threeprime[chrom][three]:
					writer.writerow([chrom, five, three, 
						threeprime[chrom][three][five]])
					num = int(threeprime[chrom][three][five])
					if num not in dist:
						dist[num] = 1
					else:
						dist[num] += 1

with open(outfilename2, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for chrom in threeprimerev:
		for three in threeprimerev[chrom]:
			if len(threeprimerev[chrom][three]) > 1:
				for five in threeprimerev[chrom][three]:
					writer.writerow([chrom, three, five, 
						threeprimerev[chrom][three][five]])
					a = int(threeprimerev[chrom][three][five])
					num = int(threeprimerev[chrom][three][five])
					if num not in dist:
						dist[num] = 1
					else:
						dist[num] += 1

print(dist)  # curious about how often junctions are being counted

# import sys, csv

# try:
# 	bed = open(sys.argv[1])
# 	outfilename = sys.argv[2]
# except:
# 	print('usage: script.py corrected_bedfile outfilename')
# 	sys.exit(1)

# threeprime = {}
# for line in bed:
# 	line = line.rstrip().split('\t')
# 	# find junctions and add to threeprime dictionary
# 	chrom = line[0]
# 	three = line[1]  
# 	five = line[2] 
# 	if chrom not in threeprime:
# 		threeprime[chrom] = {}
# 	if five not in threeprime[chrom]:
# 		threeprime[chrom][five] = {}
# 	if three not in threeprime[chrom][five]:
# 		threeprime[chrom][five][three] = 0
# 	threeprime[chrom][five][three] = line[4]

# dist = {}
# with open(outfilename, 'wt') as outfile:
# 	writer = csv.writer(outfile, delimiter='\t')
# 	for chrom in threeprime:
# 		for five in threeprime[chrom]:
# 			if len(threeprime[chrom][five]) > 1:
# 				for three in threeprime[chrom][five]:
# 					writer.writerow([chrom, three, five, threeprime[chrom][five][three]])
# 					num = int(threeprime[chrom][five][three])
# 					if num not in dist:
# 						dist[num] = 1
# 					else:
# 						dist[num] += 1
# print(dist)
