import sys#, csv

try:
	fastq = open(sys.argv[1])
	# outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py fastqfile > fastadest\n')
	sys.exit(1)

# with open(outfilename, 'wt') as outfile:
	# writer = csv.writer(outfile, delimiter='\t')
n = 0
for line in fastq:
	if n % 4 == 0:
		print('>' + line.split()[0][1:])
	elif n % 4 == 1:
		print(line.rstrip())
	n += 1