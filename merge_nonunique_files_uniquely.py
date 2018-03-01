import sys

junctions = {}
overlap = 0
for line in open(sys.argv[1]):
	line = line.rstrip().split('\t')
	j = '\t'.join(line[:3])
	junctions[j] = {}
	junctions[j]['count'] = int(line[-1])
	junctions[j]['fields'] = line[:-1]
for line in open(sys.argv[2]):
	line = line.rstrip().split('\t')
	j = '\t'.join(line[:3])
	if j in junctions:
		print('\t'.join(line[:-1]) +'\t'+str(junctions[j]['count'] + int(line[-1])))
		junctions.pop(j)
		overlap += 1
	else:
		print('\t'.join(line))
for j in junctions:
	print('\t'.join(junctions[j]['fields'])+'\t'+str(junctions[j]['count']))
sys.stderr.write(str(overlap)+'\n')