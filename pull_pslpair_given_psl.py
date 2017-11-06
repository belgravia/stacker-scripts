import sys, csv

try:
    pslquery = open(sys.argv[1], 'r')
    pslref = open(sys.argv[2], 'r')
except:
    print('pls specify 2 arguments, a psl with IDs to query in a a pair fusion psl file')
    sys.exit(1)

query = []
for line in pslquery:
    line = line.rstrip().split('\t')
    query += [line[9]]

with open(sys.argv[1][:-4] + '.pair.psl', 'wt') as pslout:
    writer = csv.writer(pslout, delimiter='\t')
    for line in pslref:
        line = line.rstrip().split('\t')
        if line[9] in query:
            writer.writerow(line)

