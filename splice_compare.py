import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	prefix = sys.argv[3]
except:
	sys.stderr.write('usage: script.py pslfile gtffile outprefix\n')
	sys.exit(1)

gtf_info = {}
transcript = False
gtf_info_3 = {}
for line in gtf:   # i really should write a gtf parser..
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand, gene = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8].split(';')[0]
	gene_id = gene[len('gene_id')+2:-1]
	if ty == 'transcript':
		transcript = True
		prev_exon_end = end
		continue
	elif ty != 'exon':
		transcript = False
		continue
	elif not transcript:
		continue

	if chrom not in gtf_info:
		gtf_info[chrom] = {}
	if prev_exon_end not in gtf_info[chrom]:
		gtf_info[chrom][prev_exon_end] = set()
	gtf_info[chrom][prev_exon_end].add(start)

	if chrom not in gtf_info_3:
		gtf_info_3[chrom] = {}
	if start not in gtf_info_3[chrom]:
		gtf_info[chrom][start] = set()
	gtf_info[chrom][start].add(prev_exon_end)

	prev_exon_end = end


fiveprime = dict.fromkeys(range(-100, 100), 0)  # 5 prime wiggle distances
threeprime = dict.fromkeys(range(-100, 100), 0)
intron = dict.fromkeys(range(-100, 100), 0)
lost5 = lost3 = 0
for line in psl:
	line = line.rstrip().split('\t')
	chrom = line[13]
	blockstarts = [int(x) for x in line[20].split(',')[:-1]]
	blocksizes = [int(x) for x in line[18].split(',')[:-1]]
	blockend = blockstarts[0] + blocksizes[0]
	for i in range(1, len(blockstarts)):
		start = blockstarts[i]
		if start - blockend >= 60:
			# assess validity of the junction
			upper = 100
			wiggle5 = wiggle3 = 0  # so for now i am not taking strand into account
			if chrom not in gtf_info:
				continue
			while blockend + wiggle5 not in gtf_info[chrom]:
				if wiggle5 == upper:
					lost5 += 1
					break
				if wiggle5 == 0:
					wiggle5 += 1
				elif wiggle5 >= 0:
					wiggle5 = wiggle5 * -1
				else:
					wiggle5 = (wiggle5-1) * -1

			if wiggle5 == upper:
				while blockend + wiggle3 not in gtf_info_3[chrom]:
					if wiggle3 == upper:
						lost3 += 1
						break
					if wiggle3 == 0:
						wiggle3 += 1
					elif wiggle3 >= 0:
						wiggle3 = wiggle3 * -1
					else:
						wiggle3 = (wiggle3-1) * -1
				if wiggle3 != upper:
					sys.stderr.write('Found one end but not the other\n')
				continue

			while start + wiggle3 not in gtf_info[chrom][blockend + wiggle5]:
				if wiggle3 == upper:
					lost3 += 1
					break
				if wiggle3 == 0:
					wiggle3 += 1
				elif wiggle3 >= 0:
					wiggle3 = wiggle3 * -1
				else:
					wiggle3 = (wiggle3-1) * -1

			if wiggle3 == upper:
				# sys.stderr.write('did not find the three prime end\n')
				continue


			wiggle = wiggle3 if abs(wiggle3) > abs(wiggle5) else wiggle5
			intron[wiggle] += 1
			# threeprime[wiggle3] += 1
			# fiveprime[wiggle5] += 1
		blockend = start + blocksizes[i]

with open(prefix + '.intron.wiggle', 'wt') as outfile_intron, open(prefix + '.5prime.wiggle', 'wt') as outfile_5, \
	open(prefix + '.3prime.wiggle', 'wt') as outfile_3:
	writer_i = csv.writer(outfile_intron, delimiter='\t')
	writer_5 = csv.writer(outfile_5, delimiter='\t')
	writer_3 = csv.writer(outfile_3, delimiter='\t')
	for i in range(-100, 100):
		writer_i.writerow([i, intron[i]])
		writer_5.writerow([i, fiveprime[i]])
		writer_3.writerow([i, threeprime[i]])
sys.stderr.write(' '.join([str(lost3), str(lost5), '\n']))