import sys
counter = 0
for line in open(sys.argv[1]):
	if counter % 4 == 0:
		print('>'+line.rstrip())
	elif counter % 4 == 1:
		print(line.rstrip())
	counter += 1
