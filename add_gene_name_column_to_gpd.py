#hello
import sys, csv
camerongtf = open('/pod/pstore/groups/brookslab/csoulette/annotations/gencode.v24.annotation.gtf', 'r')

transcript_id2gene_name = {}
for line in camerongtf:
	if line[0] == '#':
		continue
	line = line.rstrip().split('\t')
	# print(line)
	# sys.exit()
	infoline = line[8].split(';')
	if 'transcript_id' not in infoline[1]:
		continue
	transcript_id = infoline[1][infoline[1].find('ENST'):-1]
	# if transcript_id not in transcript_id2gene_name:
	gene_name = infoline[4][infoline[4].find('"')+1:-1]
	transcript_id2gene_name[transcript_id] = gene_name

filetoconvert = open(sys.argv[1], 'r')
with open(sys.argv[1][:-4] + '.genename.gpd', 'wt') as csvfile:
	writer = csv.writer(csvfile, delimiter='\t')
	i=0
	for line in filetoconvert:
		i += 1
		line = line.rstrip().split('\t')
		if 'ENST' not in line[0]:
			print(i, line)
			continue
		line = [transcript_id2gene_name[line[0]]] + line
		writer.writerow(line)
