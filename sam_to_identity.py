#!/usr/bin/env python3
# Roger Volden

'''
Usage:
python3 sam_to_identity.py alignment_with_MD.sam >out
Takes a SAM file with MD tags and returns percent identity for each read.

Output format:
readName    identity (decimal)  #matches    #mismatch   #indel
'''

import sys
import re, csv

tlengths = open(sys.argv[2])
transcript_counts = {}
sirv_sizes = {}
for line in tlengths:
    line = line.rstrip().split('\t')
    sirv_sizes[line[0]] = int(line[1])

read_transcript = {}

def readSAM(inFile):
    unaligned = 0
    for line in inFile:
        if line.startswith('@'):
            continue
        line = line.rstrip().split('\t')
        # name : matches / (matches + mismatches + insertions + deletions)
        if line[5] == '*':
            unaligned += 1
            continue
        matches = parseMD(line[-1])
        denominator, M, indel = parseCIGAR(line[5])
        mismatch = M-matches
        qname, flag, tname, pos, cigar, seq, qual = line[0], int(line[1]), line[2], int(line[3]), line[5], line[9], line[10]
        if tname not in sirv_sizes:
            continue
        identity = matches/denominator
        if qname not in read_transcript:
            read_transcript[qname] = []
        read_transcript[qname] += [(tname, matches/sirv_sizes[tname])]
        # print(line[0] + '\t' + str(identity) + '\t' + str(matches) + '\t' + str(mismatch) + '\t' + str(indel) + '\t' + target)

def parseCIGAR(cstr):
    p = re.compile(r'([MIDNSHP=X])')
    splitCstr = [i+j for i,j in zip(p.split(cstr)[::2], p.split(cstr)[1::2])]
    total, M, indel = 0, 0, 0
    for i in range(len(splitCstr)):
        if splitCstr[i][-1] in 'MID':
            total += int(splitCstr[i][:-1])
        if splitCstr[i][-1] == 'M':
            M += int(splitCstr[i][:-1])
        if splitCstr[i][-1] in 'ID':
            indel += int(splitCstr[i][:-1])
    return total, M, indel

def parseMD(mdstr):
    p = re.compile(r'[actgACTG\^]')
    return sum(int(x) for x in list(filter(None, p.split(mdstr[5:]))))

def main():
    bad = 0
    inFile = sys.argv[1]
    readSAM(open(inFile))
    num = 0
    for read in read_transcript:
        # print(read, read_transcript[read])
        num += 1
        total = sum([t[1] for t in read_transcript[read]])
        if total == 0:
            bad += 1
            continue
        for i in range(len(read_transcript[read])):
            tname = read_transcript[read][i][0]
            if tname not in transcript_counts:
                transcript_counts[tname] = 0
            transcript_counts[tname] += read_transcript[read][i][1]/total


    with open(sys.argv[3], 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for t in transcript_counts:
            writer.writerow([t, transcript_counts[t]])
    print(bad, num)
main()
