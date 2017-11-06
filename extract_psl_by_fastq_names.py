import sys, csv

try:
	psl = open(sys.argv[1])
	fastq = open(sys.argv[2])
	outfilename = sys.argv[3]
	outfilename2 = sys.argv[4]
except:
	sys.stderr.write('usage: script.py pslfile fastqfile outfilename outfilename2\n')
	sys.exit(1)

num = 0
names = set()
for line in fastq:
	if num % 4 == 0:
			names.add(line.split()[0][1:])

numlines = [0, 0]
with open(outfilename, 'wt') as outfile, open(outfilename2, 'wt') as outfile2:
	writer = csv.writer(outfile, delimiter='\t')
	writer2 = csv.writer(outfile2, delimiter='\t')
	for line in psl:
		line = line.rstrip().split('\t')
		if line[9][:line[9].find('D_')+1] in names:
			writer.writerow(line)
			numlines[0] += 1
		else:
			writer2.writerow(line)
			numlines[1] += 1
sys.stderr.write('Two files made with {} lines\n'.format(numlines))
