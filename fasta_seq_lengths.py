import sys, csv

try:
	fasta = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py fasta outfilename\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in fasta:
		line = line.rstrip()
		if line.startswith('>'):
			name = line[1:]
			continue
		writer.writerow([name, len(line)])