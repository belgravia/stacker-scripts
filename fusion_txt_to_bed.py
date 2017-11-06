import sys, csv
if len(sys.argv) > 1:
	txtfile = sys.argv[1]
else:
	txtfile = '/pod/home/alison/fusion_cll_filtered_cw_jn.txt'

if len(sys.argv) > 2:
	bedfile = sys.argv[2]
else:
	bedfile = 'fusion_cll_filtered_cw_jn.bed'

score = 0  # dummy score
with open(bedfile, 'wt') as bed:
	writer = csv.writer(bed, delimiter='\t')
	txt = open(txtfile, 'r')
	first = True
	for line in txt:
		if first:
			first = False
			#continue
		line = line.rstrip().split('\t')
		nameA, nameB = [item for item in line[1].split('--') if item]
		chrA, posA, strandA = line[2].split(':')
		chrB, posB, strandB = line[3].split(':')
		entryA = [chrA, posA, str(int(posA) + 1), nameA, score, strandA]  # dummy end position
		entryB = [chrB, posB, str(int(posB) + 1), nameB, score, strandB]
		writer.writerow(entryA)
		writer.writerow(entryB)

