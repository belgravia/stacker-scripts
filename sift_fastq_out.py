import sys, csv

try:
	fastq1 = open(sys.argv[1])
	fastq2 = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py fast1 fastq2 outfilename\n')
	sys.stderr.write('removes all instances of fastq2 names from fastq1\n')
	sys.exit(1)

n = 0
names = set()
for line in fastq2:
	if n % 4 == 0:
		names.add(line.rstrip().split()[0])
	n += 1

n = 0
removed = set()
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	name, seq, plus, qual = '', '','', ''
	for line in fastq1:
		line = line.rstrip()
		if n % 4 == 0:
			if name and name in names:
				removed.add(name)
			elif name:
				writer.writerow([name])
				writer.writerow([seq])
				writer.writerow([plus])
				writer.writerow([qual])
			name = line.split()[0]
		elif n % 3 == 0:
			qual = line
		elif n % 2 == 0:
			plus = line
		else:
			seq = line
		n += 1
	if not name in names:
		writer.writerow([name])
		writer.writerow([seq])
		writer.writerow([plus])
		writer.writerow([qual])

if removed and len(removed) < 10:
	print(removed)
else:
	sys.stderr.write('{} entries removed\n'.format(len(removed)))
