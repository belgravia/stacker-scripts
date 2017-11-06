import sys, csv 
import numpy

try:
	fastq = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py fastq outfilename [read]\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	n = 0
	for line in fastq:
		if n % 4 == 0:
			name = line.rstrip().split()[0][1:]
		if n % 4 == 3:
			line = line.rstrip()
			scores = []
			for char in line:
				scores += [ord(char)]
			avg = sum(scores)/float(len(scores))
			sd = numpy.std(scores)
			writer.writerow([name, avg, sd])
		n += 1