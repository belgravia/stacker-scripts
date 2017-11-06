import sys, csv

try:
	psl = open(sys.argv[1])
	# outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py bedfile > newbedfile\n')
	sys.exit(1)

prevline = ''
for line in psl:
	line = line.rstrip().split('\t')
	if not prevline:
		prevline = line
		continue
	if line[0] + line[1] + line[2] != prevline[0] + prevline[1] + prevline[2]:
		print('\t'.join(prevline))
		prevline = line
	else:
		line[4] = str(int(line[4]) + int(prevline[4]))
print('\t'.join(line))
