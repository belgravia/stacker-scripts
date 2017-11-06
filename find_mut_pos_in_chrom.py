import sys, csv

try:
	gtf = open(sys.argv[1])
	gene = sys.argv[2]
	aa_pos = int(sys.argv[3]) * 3
except:
	sys.stderr.write('usage: script.py gtf gene aa_pos\n')
	sys.exit(1)

transcript = ''
cds_len = 0
for line in gtf:
	if line.startswith('#'):
		continue
	line = line.rstrip().split('\t')
	info = line[8]

	if gene in info and line[2] == 'CDS':
		this_transcript = info[info.find('ENST'):]
		this_transcript = this_transcript[:this_transcript.find('"')]
		cds_start, cds_end = int(line[3]), int(line[4])
		if this_transcript != transcript:
			if cds_len != 0:
				print('{} total CDS length: {}'.format(transcript, cds_len))
			transcript = this_transcript
			cds_len = cds_end - cds_start
		else:
			if cds_len + cds_end - cds_start >= aa_pos and cds_len < aa_pos:
				pos = cds_end - (aa_pos - cds_len)   # it's on the minus strand. 
				# pos = cds_start + aa_pos - cds_len   # use this if the gene is on the pos strand
				print('{} {}'.format(transcript, pos))
			cds_len += cds_end - cds_start


