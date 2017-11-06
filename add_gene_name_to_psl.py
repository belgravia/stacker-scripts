#hello
import sys, csv

try:
	psl = open(sys.argv[1]) # my step 1 potential fusion psl
except:
	print('1 required argument, a psl')
	sys.exit(1)
try:
	if len(sys.argv) > 2:  # optionally specify a gtf name
		smallergtf = open(sys.argv[2], 'r')
	else:
		smallergtf = open('gencodev25_cathyfusionnames_intersect.gtf', 'r')
except:
	print('bad gtf path')
	sys.exit(1)

if len(sys.argv) > 3:  # optionally specify you want to keep unmapped psl entries too
	if 'add' in sys.argv[3]:
		add = True
		minlen = int(sys.argv[3][3:])
		if not minlen:
			minlen = 100
	else:
		minlen = int(sys.argv[3])
		add = False
else:
	add = False
	minlen = 100

def entry_gene_map(entry, line):
	if entry[13] != line[0]:  # compare chromosomes
		return False
	# compare within range of gene
	if (10 + int(line[3])) <= int(entry[15]) and int(entry[16]) <= (int(line[4]) + 10):
		return True
	return False

gtfentries = []
for line in smallergtf:
	if line[0] == '#':
		continue
	line = line.rstrip().split('\t')
	if line[2] != 'gene':
		continue
	gtfentries += [line]
setofchr = set([entry[0] for entry in gtfentries])
entries_to_keep = []
printonce = True

for entry in psl:
	entry = entry.rstrip().split('\t')
	if entry[13] not in setofchr or int(entry[0]) < minlen:
		continue
	added = False
	# print('asdf')
	for line in gtfentries:
		infoline = line[8].split(';')
		gene_name = infoline[3][infoline[3].find('"')+1:-1]
		if entry_gene_map(entry, line):  # psl entry and gtf line
			if gene_name not in entry[-1]:
				entry[-1] += gene_name + ';'
				added = True
				if printonce:
					print('it works :D')
					printonce = False
	if added or add:
		entries_to_keep += [entry]
		

# write modified psl entries to a file
# probably should write on the fly
with open(sys.argv[1][:-4] + '.genename.psl', 'wt') as csvfile:
	writer = csv.writer(csvfile, delimiter='\t')
	for entry in entries_to_keep:
		writer.writerow(entry)

with open('done', 'wt') as csvfile:
	writer = csv.writer(csvfile, delimiter='\t')
	writer.writerow(['hi'])


