import sys, csv

try:
	psl = open(sys.argv[1])
	if len(sys.argv) > 2:
		sirvs = True
	else:
		sirvs = False
except:
	sys.stderr.write('usage: script.py pslfile [sirvs?] > outfilename\n')
	sys.exit(1)

for line in psl:
	line = line.rstrip().split('\t')
	chrom, strand, score, name = line[13], line[8], line[0], line[9]
	tstarts = [int(n) for n in line[20].split(',')[:-1]]  # target starts
	bsizes = [int(n) for n in line[18].split(',')[:-1]]  # block sizes
	
	# gene_id = name[:name.rfind('-')] mandalorion specific
	#transcript_id = gene_id + '_' + name[name.find('Isoform'):]
	if 'ENSG' in name:
		gene_id = name[name.find('ENSG'):]
	else:
		gene_id = name[name.find('chr'):]
	gene_id = gene_id[:gene_id.find('_')]
	transcript_id = name[:name.rfind('_')]	

	if sirvs:
		gene_id = chrom
		tid = name[name.find('-')+1:name.find('_')]
		if len(str(tid)) == 1:
			transcript_id = gene_id + '00' + str(tid)
		elif len(str(tid)) == 2:
			transcript_id = gene_id + '0' + str(tid)
		elif len(str(tid)) == 3:
			transcript_id = gene_id + str(tid)
		else:
			sys.stderr.write('too many isoforms found for gene (1000+)\n')
			sys.exit()
	for b in range(len(tstarts)):
		exon_assignment = transcript_id + '_' + str(b)
		endstring = 'gene_id \"{}\"; transcript_id \"{}\"; exon_assignment \"{}\";'\
				.format(gene_id, transcript_id, exon_assignment)
		print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[b]), \
			str(tstarts[b]+bsizes[b]), '.', strand, str(score), endstring]))
