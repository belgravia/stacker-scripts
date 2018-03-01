import sys, csv

try:
	fastq = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py fastqfile outfilename\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	n = 0
	removed = -1
	name, seq, plus, qual = '', '','', ''
	for line in fastq:
		line = line.rstrip()
		if n % 4 == 0:
			if not name or len(seq)!=len(qual):
				removed += 1
			else:
				writer.writerow([name])
				writer.writerow([seq])
				writer.writerow([plus])
				writer.writerow([qual])
			name = line
		elif n % 3 == 0:
			qual = line
		elif n % 2 == 0:
			plus = line
		else:
			seq = line
		n += 1
	if name and len(seq) == len(qual):
		writer.writerow([name])
		writer.writerow([seq])
		writer.writerow([plus])
		writer.writerow([qual])
	else:
		removed += 1

sys.stderr.write('{} entries removed\n'.format(removed))
