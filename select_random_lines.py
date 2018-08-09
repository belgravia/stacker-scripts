import sys, random

try:
	fle = open(sys.argv[1])
	percent = float(sys.argv[2])
except:
	# sys.stderr.write('usage: script.py filename #lines > subset\n')
	sys.stderr.write('usage: script.py filename percent_lines > subset\n')
	sys.exit()

for line in fle:
	# if percent < 1:
	if random.random() < percent:
		print(line.rstrip())
