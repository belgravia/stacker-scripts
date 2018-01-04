import sys, csv, re

try:
	sam = open(sys.argv[1])
	# tlengths = open(sys.argv[2])
	outfilename = sys.argv[2]
except:
	sys.stderr.write('usage: script.py samfile tableout\n')
	sys.exit(1)

transcript_lengths = {}
read_transcript = {}
transcript_counts = {}

# for line in tlengths:
# 	line = line.rstrip().split('\t')
# 	transcript_lengths[line[0]] = int(line[1])

for line in sam:
	if line.startswith('@'):
		continue
	line = line.rstrip().split('\t')
	qname, flag, tname, pos, cigar, seq, qual = line[0], int(line[1]), line[2], int(line[3]), line[5], line[9], line[10]
	if qname not in read_transcript:
		read_transcript[qname] = []
	if tname == '*':
		continue
	read_transcript[qname] += [tname]

for read in read_transcript:
	for transcript in read_transcript[read]:
		if transcript not in transcript_counts:
			transcript_counts[transcript] = 0
		transcript_counts[transcript] += 1/float(len(read_transcript[read]))

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for transcript in transcript_counts:
		if transcript == '*':
			continue
		writer.writerow([transcript, transcript_counts[transcript]])
