import sys, csv

try:
	grep = sys.argv[1]
	fle = open(sys.argv[2])
except:
	sys.stderr.write('usage: script.py grepstring file\n')
	sys.exit(1)

num = 0
for line in fle:
	if grep in line:
		print('{}: {}'.format(num, line.rstrip()))
	num += 1
