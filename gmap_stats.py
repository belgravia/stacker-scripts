import sys, csv

try:
	fastq = open(sys.argv[1])
	psl = open(sys.argv[2])
except:
	sys.stderr.write('usage: script.py fastq pslfile\n')
	sys.exit(1)

readnames = {}
i = 0
for line in fastq:
	if i % 4 == 0:
		readnames[line.split()[0][1:]] = 0
	i += 1


mm = 0  # multiply mapping read counts
mapped = 0  # at least one alignment 
for line in psl:
	line = line.rstrip().split('\t')
	line[9] = line[9]#[:-4]
	readnames[line[9]] += 1
	if readnames[line[9]] == 2:
		mm += 1
	elif readnames[line[9]] == 1:
		mapped += 1

sys.stderr.write('mapped/total = {}/{} = {} % aligned\n'.format(mapped, len(readnames), 100*float(mapped)/len(readnames)))
sys.stderr.write('{} multiply mapping reads\n'.format(mm))
