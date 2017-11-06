import sys, csv

try:
	gene_list = open(sys.argv[1])
	gtf_in = open(sys.argv[2])
	gtf_out = sys.argv[3]
except:
	sys.stderr.write('usage: script.py gene_list gtf_in gtf_out\n')
	sys.exit(1)

gl = set()
for line in gene_list:
	gl.add(line.rstrip())

with open(gtf_out, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in gtf_in:   # i really should write a gtf parser..
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		chrom, ty, start, end, strand, gene = line[0], line[2], int(line[3]), int(line[4]), line[6], line[8].split(';')[0]
		gene_id = gene[len('gene_id')+2:-1]
		if gene_id in gl:
			writer.writerow(line)
