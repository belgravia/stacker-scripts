import sys, re

try:
	sam = open(sys.argv[1])
except:
	sys.stderr.write('usage: script.py samfile > outfilename\n')
	sys.exit(1)

for line in sam:
	if line.startswith('@'):
		print(line.rstrip())
		continue
	line = line.rstrip().split('\t')
	seq, qual, cigar = line[9], line[10], line[5]
	matches = re.findall('([0-9]+)([A-Z])', cigar)
	matchlen = 0
	for m in matches:
		if m[1] in ['M', 'D']:  # both matches and deletions consume reference
			matchlen += int(m[0])
	if cigar != '*':
		line[5] = str(matchlen) + 'M'
		if len(seq) > matchlen:
			seq = seq[:matchlen]
			qual = qual[:matchlen]
		elif len(seq) < matchlen:
			seq += 'N'*(matchlen - len(seq))
			qual += '3'*(matchlen - len(qual))
		line[9] = seq
		line[10] = qual
	print('\t'.join(line))

