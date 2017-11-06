import sys, csv

try:
	psl = open(sys.argv[1])
	col = int(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py pslfile column_number outfilename\n')
	sys.stderr.write('for a format suitable for goseq r package\n')
	sys.exit(1)

genes = {}
for line in psl:
	line = line.rstrip().split('\t')
	info = line[col]
	if line[1] != '1':  # task-specific
		continue
	if float(line[6]) > 0:
		wtgreater = 1
	else:
		wtgreater = 0
	while 'ENSG' in info:
		gene = info[info.find('ENSG'):info.find('.')]
		if not gene in genes or genes[gene]['ir'] == '0':
			if gene not in genes:
				genes[gene] = {}
			genes[gene]['ir'] = line[1]
			genes[gene]['wt'] = wtgreater
		info = info[info.find('-')+1:]


with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for g in genes:
		writer.writerow([g, genes[g]['wt']])

