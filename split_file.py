import sys, csv

try:
	psl = open(sys.argv[1])
	outfilenamebase = sys.argv[2]
	lines = int(sys.argv[3])
except:
	sys.stderr.write('usage: script.py pslfile outfilenamebase numlines\n')
	sys.exit(1)

count = 0
base = 0
line = psl.readline().rstrip().split('\t')
while line != ['']:
	outfilename = outfilenamebase + '_' + str(base) + '.psl'
	with open(outfilename, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t')
		while count != lines:
			writer.writerow(line)
			count += 1
			line = psl.readline().rstrip().split('\t')
			if line == ['']:
				break
	base += 1
	count = 0