import sys, random

try:
	fle = open(sys.argv[1])  # psl
	percent = float(sys.argv[2])  # what % of lines to randomly keep
except:
	sys.stderr.write('usage: script.py filename percent_lines > subset\n')
	sys.exit()

for line in fle:
	if random.random() < percent:
		print(line.rstrip())
