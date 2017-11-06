import sys

sys.stderr.write('usage: script.py [some number of files] > catted\n')

for fle in sys.argv[1:]:
	for line in open(fle):
		print(line.rstrip())

