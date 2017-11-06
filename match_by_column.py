import sys, csv

try:
	files = []
	colmatch = []
	colpull = []
	for i in range(1, len(sys.argv)-1, 3):
		files += [sys.argv[i]]
		colmatch += [int(sys.argv[i+1])]
		colpull += [int(sys.argv[i+2])]
	if len(sys.argv) % 3 != 2:
		1/0
	outfilename = sys.argv[-1]
except:
	sys.stderr.write('usage: script.py [filename colmatch colpull] [...] outfilename\n')
	sys.exit(1)

alldata = {}
for i in range(len(files)):
	f = files[i]
	colm = colmatch[i]
	colp = colpull[i]
	for line in open(f):
		line = line.rstrip().split('\t')
		if line[colm] in alldata:
			if len(alldata[line[colm]])+1 == i:
				continue
			else:
				alldata[line[colm]] += ['NA']*(len(alldata[line[colm]])-i)
				alldata[line[colm]] += [line[colp]]
		else:
			alldata[line[colm]] = ['NA']*(i) + [line[colp]]

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for d in alldata:
		if len(alldata[d]) != len(files):
			alldata[d] += ['NA']*(len(files)-len(alldata[d]))
		writer.writerow([d] + alldata[d])

