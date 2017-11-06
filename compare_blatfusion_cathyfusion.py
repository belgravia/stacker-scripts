from __future__ import print_function
import sys, csv
def myRound(num, sigfigs=2, maxnumzeros=7):
	numdigits = len(num) if len(num) > maxnumzeros else maxnumzeros
	roundfactor = 10 ** (numdigits - sigfigs)
	return str(int(num)//roundfactor * roundfactor)
	# return '0'

def validFusion(s1, e1, s2, e2):
	s1, e1 = sorted([s1, e1], key=int)
	s2, e2 = sorted([s2, e2], key=int)
	s1, e1, s2, e2 = [int(i) for i in [s1, e1, s2, e2]]
	if s1 - 10 < s2 and s1 + 10 > s2:  # if s2 is +/- 10 of s1
		return False
	if (s2 < e1 - 10 and s2 > s1) or (s1 < e2 - 10 and s1 > s2):  #  any more overlap than 10 is bad
		return False
	return True
#7b0b5eea-6a9c-4a39-9176-9407f5012d48_Basecall_2D_2d
# one required argument

def compare_entries_fusion_not_stringent(chrA, chrB):
	for fus in cathy_fusions:
		if fus == (chrA + '--' + chrB) or fus == (chrB + '--' + chrA):
			return True
	return False

def compare_entries_fusion(posA, posB):
	for fus in cathy_fusions:
		# if posA == fus
		continue

def process_entries(entries, line, name):
		# entries.sort(key=lambda x: x[13])  # target sequence name is the tie breaker
		# entries.sort(key=lambda x: int(x[0]), reverse=True)  # number matches
		# entries = [e for e in entries if int(e[0]) > 50]
		if not name or len(entries) < 2:  # init
			return
		if len(entries) > 2:
			print('panic')
		chrA = entries[0][13] + ':' + myRound(entries[0][15])
		chrB = entries[1][13] + ':' + myRound(entries[1][15])
		if myRound(entries[0][15]) != myRound(entries[0][16]) or myRound(entries[1][15]) != myRound(entries[1][16]):
			print(myRound(entries[0][15]), myRound(entries[0][16]))
			print(myRound(entries[1][15]), myRound(entries[1][16]))
		chrA, chrB = sorted([chrA, chrB])
		key = (chrA, chrB)
		if compare_entries_fusion_not_stringent(chrA, chrB):
			print('added in dict')
			if key not in fusionfreq:
				fusionentries[key] = entries
			else:
				fusionentries[key] += entries	
try:
	bed = open(sys.argv[1])
	psl = open(sys.argv[2]) # my step 1 potential fusion psl
except:
	print('2 required arguments, cathy hg38 bedfile and my blat fusion psl')
	sys.exit(1)

cathy_fusions = []
name = ''
for line in bed:
	line = line.split('\t')
	if name:
		cathy_fusions += [name + '--' + line[0] + ':' + myRound(line[1])]
		name = ''
	else:
		name = line[0] + ':' + myRound(line[1])

name = ''
entries = []
fusionfreq = {}
fusionentries = {}
for line in psl:
	line = line.rstrip().split('\t')
	if line[9] != name:
		if len(entries) > 2:
			print(entries)	
			if entries[1] == entries[3] and entries[0] == entries[2]:
				entries = entries[:2]
		process_entries(entries, line, name)
		name = line[9]
		entries = [line]
		continue
	entries += [line]
process_entries(entries, line, name)

if not fusionentries:
	print('did not work D:')
with open(sys.argv[2][:-4] + '_fusions_in_agreement.psl', 'wt') as csvfile:
	writer = csv.writer(csvfile, delimiter='\t')
	for k in fusionentries:
		for entry in fusionentries[k]:
			writer.writerow(entry)
