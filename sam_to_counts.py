import sys, csv, re

try:
	sam = open(sys.argv[1])
	tlengths = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	sys.stderr.write('usage: script.py samfile tlengths table\n')
	sys.exit(1)

transcript_lengths = {}
read_transcript = {}
transcript_counts = {}

for line in tlengths:
	line = line.rstrip().split('\t')
	transcript_lengths[line[0]] = int(line[1])

for line in sam:
	if line.startswith('@'):
		print(line)
		continue
	line = line.rstrip().split('\t')
	qname, flag, tname, pos, cigar, seq, qual = line[0], int(line[1]), line[2], int(line[3]), line[5], line[9], line[10]
	if qname not in read_transcript:
		read_transcript[qname] = []
	print(qname, tname)
	matches = re.findall('([0-9]+)([A-Z])', cigar)
	match_count = 0
	for m in matches:
		if m[1] == 'M':
			match_count += int(m[0])
	if tname == '*':
		continue
	read_transcript[qname] += [(tname, match_count/float(transcript_lengths[tname]))]

for read in read_transcript:
	print(read_transcript[read])
	total = sum([t[1] for t in read_transcript[read]])
	for i in range(len(read_transcript[read])):
		transcript = read_transcript[read][i][0]
		if transcript not in transcript_counts:
			transcript_counts[transcript] = 0

		transcript_counts[transcript] += read_transcript[read][i][1]/float(total)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	# for read in read_transcript:
	# 	if read_transcript[read] == ['*']:
	# 		continue
	# 	if len(read_transcript[read]) != len(set(read_transcript[read])):
	# 		print(read)
	# 		print(read_transcript[read])
	# 		bad  += 1
	# 	if len(read_transcript[read]) > 1:
	# 		continue
	# 	writer.writerow([read] + read_transcript[read])
	for transcript in transcript_counts:
		if transcript == '*':
			continue
		writer.writerow([transcript, transcript_counts[transcript]])
