import sys, csv

try: 
	pslfile = sys.argv[1]
except:
	print('need 2 arguments, psl and optionally the bed output file name')

if len(sys.argv) > 2:
	bedfile = sys.argv[2]
else:
	bedfile = sys.argv[1][:-3] + 'bed'

with open(bedfile, 'wt') as bed:
	writer = csv.writer(bed, delimiter='\t')
	psl = open(pslfile, 'r')
	for line in psl:
		line = line.rstrip().split('\t')
		writer.writerow([line[13], line[15], line[16], line[9], line[0], line[8]])

