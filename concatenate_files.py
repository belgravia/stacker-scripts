import sys

sys.stderr.write('usage: script.py [some number of files] > catted\n')

header = False
for fle in sys.argv[1:]:
	for line in open(fle):
		if not line.startswith('@'):
			print(line.rstrip())
		elif not header:
			print(line.rstrip())
	header = True
