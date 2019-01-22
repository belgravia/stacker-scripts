import sys

if len(sys.argv) < 2:
	sys.stderr.write('usage: script.py table1 table2 table3... > combined.txt \n')
	sys.stderr.write('matched by the first column\n')
	sys.exit()

combined = {}
header = ['gene']
for i in range(len(sys.argv[1:])):
	header += [sys.argv[i+1]]
	for line in open(sys.argv[i+1]):
		line = line.rstrip().split('\t')
		if line[0] not in combined:
			combined[line[0]] = [0]*i + [line[1]]
		else:
			combined[line[0]] += [line[1]]

print('\t'.join(header))

genes = sorted(list(combined.keys()))
for g in genes:
	print('\t'.join([g] + combined[g]))