import sys, csv

try:
	sam = open(sys.argv[1])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py samfile table\n')
	sys.exit(1)

read_transcript = {}

for line in sam:
	if line.startswith('@'):
		continue
	line = line.rstrip().split('\t')
	qname, flag, tname, pos, cigar, seq, qual = line[0], int(line[1]), line[2], int(line[3]), line[5], line[9], line[10]
	if qname not in read_transcript:
		read_transcript[qname] = []
	read_transcript[qname] += [tname]

bad=0
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for read in read_transcript:
		if read_transcript[read] == ['*']:
			continue
		if len(read_transcript[read]) != len(set(read_transcript[read])):
			print(read)
			print(read_transcript[read])
			bad  += 1
		if len(read_transcript[read]) > 1:
			continue
		writer.writerow([read] + read_transcript[read])

print(bad)