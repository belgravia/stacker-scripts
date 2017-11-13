import sys, csv

try:
	psl = open(sys.argv[1])
	gtf = open(sys.argv[2])
	prefix = sys.argv[3]
	if len(sys.argv) > 4:
		header = sys.argv[4]
	else:
		header = ''
except:
	sys.stderr.write('usage: script.py pslfile gtffile outprefix [noheader] \n')
	sys.stderr.write('first make sure the psl\'s strand is correct\n')
	sys.exit(1)

gtf_info = {'+': {}, '-': {}}
gtf_info_3 = {'+':{}, '-':{}}
prev_transcript = ''
for line in gtf:   # i really should write a gtf parser..
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	chrom, ty, start, end, strand, info = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8].split(';')
	gene_id = info[0][len('gene_id')+2:-1]
	transcript_id = info[1][len('transcript_id')+2:-1]

	if ty != 'exon':
		continue
	if transcript_id != prev_transcript:
		prev_transcript = transcript_id
		prev_exon_end = end
		continue

	start -= 1
	if chrom not in gtf_info[strand]:
		gtf_info[strand][chrom] = {}
	if prev_exon_end not in gtf_info[strand][chrom]:
		gtf_info[strand][chrom][prev_exon_end] = set()
	gtf_info[strand][chrom][prev_exon_end].add(start)

	if chrom not in gtf_info_3[strand]:
		gtf_info_3[strand][chrom] = {}
	if start not in gtf_info_3[strand][chrom]:
		gtf_info_3[strand][chrom][start] = set()
	gtf_info_3[strand][chrom][start].add(prev_exon_end)

	prev_exon_end = end

print(gtf_info)
print(gtf_info_3)

fiveprime = dict.fromkeys(range(-100, 100), 0)  # 5 prime wiggle distances
threeprime = dict.fromkeys(range(-100, 100), 0)
intron = dict.fromkeys(range(-100, 100), 0)
lost5 = lost3 = 0
if header:
	psl.readline()
for line in psl:
	line = line.rstrip().split('\t')
	chrom = line[13]
	blockstarts = [int(x) for x in line[20].split(',')[:-1]]
	blocksizes = [int(x) for x in line[18].split(',')[:-1]]
	if len(blockstarts) == 1:  # single exon transcript has no junctions
		continue
	blockend = blockstarts[0] + blocksizes[0]
	strand = line[8]
	for i in range(1, len(blockstarts)):
		start = blockstarts[i]

		if start - blockend >= 60:
			# assess validity of the junction
			upper = 100
			wiggle5 = wiggle3 = 0  # so for now i am not taking strand into account
			if chrom not in gtf_info[strand]:
				continue
			while blockend + wiggle5 not in gtf_info[strand][chrom]:  # find the nearest annotated 5'
				if wiggle5 == upper:
					lost5 += 1
					break
				if wiggle5 == 0:
					wiggle5 += 1
				elif wiggle5 >= 0:
					wiggle5 = wiggle5 * -1
				else:
					wiggle5 = (wiggle5-1) * -1

			while start + wiggle3 not in gtf_info_3[strand][chrom]:  # find the nearest annotated 3'
				if wiggle3 == upper:
					lost3 += 1
					break
				if wiggle3 == 0:
					wiggle3 += 1
				elif wiggle3 >= 0:
					wiggle3 = wiggle3 * -1
				else:
					wiggle3 = (wiggle3-1) * -1

			if wiggle5 == upper or wiggle3 == upper:
				sys.stderr.write('Junction {}{}:{}-{} does not match any annotated junctions\n'.format(chrom,
					strand, blockend, start))
				continue

			if start + wiggle3 not in gtf_info[strand][chrom][blockend + wiggle5]:
				wiggle3_new = 0  # remap, it is possible there is a bigger wiggle that will allow the 5' and 3' agree
				while start + wiggle3_new not in gtf_info[strand][chrom][blockend + wiggle5]:
					if wiggle3_new == upper:
						break
					if wiggle3_new == 0:
						wiggle3_new += 1
					elif wiggle3_new >= 0:
						wiggle3_new = wiggle3_new * -1
					else:
						wiggle3_new = (wiggle3_new-1) * -1

				if wiggle3_new == upper:
					wiggle5 = 0
					while blockend + wiggle5 not in gtf_info_3[strand][chrom][start+wiggle3]:  # find the nearest annotated 5'
						if wiggle5 == upper:
							break
						if wiggle5 == 0:
							wiggle5 += 1
						elif wiggle5 >= 0:
							wiggle5 = wiggle5 * -1
						else:
							wiggle5 = (wiggle5-1) * -1
					if wiggle5 == upper:
						sys.stderr.write("Junction {}{}:{}-{} matches to different annotated junctions\n".format(chrom,
							strand, blockend, start))
						continue
				else:
					wiggle3 = wiggle3_new

			if strand == '-':
				wiggle5, wiggle3 = -1*wiggle3, -1*wiggle5

			wiggle = wiggle3 if abs(wiggle3) > abs(wiggle5) else wiggle5
			intron[wiggle] += 1
			threeprime[wiggle3] += 1
			fiveprime[wiggle5] += 1
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