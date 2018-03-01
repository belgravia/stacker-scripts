import sys, csv

try:
	fastq = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py fastq outfilename\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	n = 0
	removed = 0
	name, seq, plus, qual = '', '','', ''
	for line in fastq:
		line = line.rstrip().split('\t')
		if n % 4 == 0:
			if name and (not seq or not qual):
				removed += 1
			else:
				writer.writerow(name)
				writer.writerow(seq)
				writer.writerow(plus)
				writer.writerow(qual)
			name = line
		elif n % 3 == 0:
			qual = line
		elif n % 2 == 0:
			plus = line
		else:
			seq = line
		n += 1

sys.stderr.write('removed {} entries\n'.format(removed))
