from __future__ import print_function
import sys, csv
def myRound(num, sigfigs=3, maxnumzeros=7):
	numdigits = len(num) if len(num) > maxnumzeros else maxnumzeros
	roundfactor = 10 ** (numdigits - sigfigs)
	return str(int(num)//roundfactor * roundfactor)

def compare_entries_fusion_not_stringent(chrA, chrA2, chrB, chrB2):
	for fus in idp_fusions:
		fusA, empty, fusB = fus.split('-')
		if (fusA in [chrA, chrA2] and fusB in [chrB, chrB2]) or \
		(fusB in [chrA, chrA2] and fusA in [chrB, chrB2]):
			return True
	return False

def compare_entries_fusion(posA, posB):
	for fus in idp_fusions:
		# if posA == fus
		continue

def process_entries(entries, line, name):
		# entries.sort(key=lambda x: x[13])  # target sequence name is the tie breaker
		# entries.sort(key=lambda x: int(x[0]), reverse=True)  # number matches
		# entries = [e for e in entries if int(e[0]) > 50]
		if not name or len(entries) < 2:  # init
			return
		chrA = entries[0][13] + ':' + myRound(entries[0][15])  # e.g. chr1:103000
		chrA2 = entries[0][13] + ':' + myRound(entries[0][16])

		chrB = entries[1][13] + ':' + myRound(entries[1][15])
		chrB2 = entries[1][13] + ':' + myRound(entries[1][16])

		# if myRound(entries[0][15]) != myRound(entries[0][16]):
			# print(myRound(entries[0][15]), myRound(entries[0][16]))
		# if myRound(entries[1][15]) != myRound(entries[1][16]):
			# print(myRound(entries[1][15]), myRound(entries[1][16]))
		chrA_, chrB_ = sorted([chrA, chrB])
		key = (chrA_, chrB_)
		if compare_entries_fusion_not_stringent(chrA, chrA2, chrB, chrB2):
			print('added in dict')
			if key not in fusionfreq:
				fusionentries[key] = entries
			else:
				fusionentries[key] += entries	
try:
	psl = open(sys.argv[1]) # my step 1 potential fusion psl
	tsv = open(sys.argv[2])
except:
	print('2 required arguments, my blat fusion psl and idp output tsv')
	sys.exit(1)

try:
	infobool = sys.argv[3]
	infobool = True
except:
	infobool = False

idp_fusions = []
if not infobool:
	for line in tsv:
		line = line.split('\t')
		idp_fusions += [line[1]+':'+myRound(line[2])+'--'+line[5]+':'+myRound(line[6])]
else:
	for line in tsv:  # info file
		line = line.rstrip().split('\t')
		left, right = line[1].split('/')
		left = left.split('_')
		right = right.split('_')
		try:
			idp_fusions += [left[0]+':'+myRound(left[1])+'--'+right[0]+':'+myRound(right[1])]
		except:
			continue

print(idp_fusions)

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
with open(sys.argv[1][:-4] + '.idp.psl', 'wt') as csvfile:
	writer = csv.writer(csvfile, delimiter='\t')
	for k in fusionentries:
		for entry in fusionentries[k]:
			writer.writerow(entry)
