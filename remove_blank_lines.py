import sys, csv

try:
	psl = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py pslfile outfilename\n')
	sys.exit(1)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in psl:
		line = line.rstrip()
		if line:
			writer.writerow([line])