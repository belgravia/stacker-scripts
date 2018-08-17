import sys, csv

# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.stats as stats

try:
    file1 = open(sys.argv[1])
    outfilename = sys.argv[2]
    somebool = True if len(sys.argv)>3 else False
except:
    print('usage: script.py')
    sys.exit()

def get_junctions_psl(starts, sizes):
    junctions = set()
    if len(starts) != 1:
        for b in range(len(starts)-1):
            junctions.add((starts[b]+sizes[b], starts[b+1]))
        return junctions

for line in psl:
    # line = line.rstrip().split()

    # if not line.startswith('@'):  # sam
    #   qname, flag, tname, pos, cigar, seq, qual = line[0], int(line[1]), line[2], int(line[3]), line[5], line[9], line[10]
    
    # if not line.startswith('#'):  # gtf
    #   chrom, ty, start, end, strand = line[0], line[2], int(line[3]), int(line[4]), line[6]
    #   # transcript_id = line[8][line[8].find('transcript_id')+15:line[8].find('gene_type')-3]  # p specific to gencode v24
    #   transcript_id = line[8][line[8].find('transcript_id')+15:]
    #   transcript_id = transcript_id[:transcript_id.find('"')]
    #   if ty != 'exon':
    #       continue

    chrom, name, start, end = line[13], line[9], int(line[15]), int(line[16])  # psl
    sizes = [int(n) for n in line[18].split(',')[:-1]]
    starts = [int(n) for n in line[20].split(',')[:-1]]



with open(outfilename, 'wt') as outfile:
    writer = csv.writer(outfile, delimiter='\t')