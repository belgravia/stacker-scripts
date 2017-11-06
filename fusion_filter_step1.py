from __future__ import print_function
import sys, csv

def myRound(num):
	numdigits = len(num) if len(num) > 6 else 6
	roundfactor = 10 ** (numdigits - 4)
	return str(int(num)//roundfactor * roundfactor)

def validFusion(s1, e1, s2, e2):
	s1, e1 = sorted([s1, e1], key=int)
	s2, e2 = sorted([s2, e2], key=int)
	s1, e1, s2, e2 = [int(i) for i in [s1, e1, s2, e2]]
	if s1 - 10 < s2 and s1 + 10 > s2:  # if s2 is +/- 10 of s1
		return False
	if (s2 < e1 - 10 and s2 > s1) or (s1 < e2 - 10 and s1 > s2):  #  any more overlap than 10 is bad
		return False
	if (s2 - e1 > 500 and s2 > s1) or (s1 - e2 > 500 and s1 > s2):
		return False
	return True
#7b0b5eea-6a9c-4a39-9176-9407f5012d48_Basecall_2D_2d

def process_entries(entries, line, name):
		entries.sort(key=lambda x: x[13])  # target sequence name is the tie breaker
		entries.sort(key=lambda x: int(x[0]), reverse=True)  # number matches
		entries = [e for e in entries if int(e[0]) > 150]
		if not name or len(entries) < 2:
			return
		for pf1 in entries[:2]:  # potential fusion part 1
			for pf2 in entries[1:]:
				if validFusion(pf1[11], pf1[12], pf2[11], pf2[12]): 
					chra_start = pf1[13] + pf1[8] + ':' + myRound(pf1[15])
					chra_end = pf1[13] + pf1[8] + ':' + myRound(pf1[16])
					chrb_start = pf2[13] + pf2[8] + ':' + myRound(pf2[15])
					chrb_end = pf2[13] + pf2[8] + ':' + myRound(pf2[16])
					# record += [pf1, pf2]
					return [pf1, pf2]
print('hi')
try:
	psl = open(sys.argv[1])
except:
	print('1 required argument, a psl file to be filtered by fusion')
	sys.exit(1)
name = ''
entries = []
record = []
for line in psl:
	line = line.rstrip().split('\t')
	if line[9] != name:
		result = process_entries(entries, line, name)
		if result:
			record += result
		name = line[9]
		entries = [line]
		continue
	entries += [line]
result = process_entries(entries, line, name)
if result:
	record += result

with open(sys.argv[1][(sys.argv[1].rfind('/')+1):-4] + '_fusions_step1.psl', 'wt') as csvfile:
	writer = csv.writer(csvfile, delimiter='\t')
	for line in record:
		writer.writerow(line)
