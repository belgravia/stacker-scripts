import sys, csv, random

# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.stats as stats

try:
	txt = open(sys.argv[1])
	outfilename = sys.argv[2]
	if len(sys.argv) > 3:
		bounds = sys.argv[3]
		incr =  int(bounds[bounds.rfind(','):])
		upper = int(bounds[bounds.find(',')+1:bounds.rfind(',')])
		lower = int(bounds[:bounds.find(',')])
	else:
		lower = 10
		upper = 100
		incr = 10
	if len(sys.argv) > 4:
		replicates = int(sys.argv[4])
	else:
		replicates = 3
except:
	sys.stderr.write('usage: script.py txtfile outtxt [bounds] [replicates]\n')
	sys.stderr.write('txtfile format: event_name\tevent_result\n')
	sys.stderr.write('bounds format in percentage: lower_result,upper_result,increment\n')
	sys.stderr.write('e.g. default bounds are 10,100,10\n')
	sys.exit()

event_results = []
f = 0
for line in txt:
	line = line.rstrip().split('\t')
	event_results += [line[1]]
	f += 1
	if f == 10000:
		break
num_events = len(event_results)

with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	writer.writerow(['percentage', 'result'])
	for percent in range(lower, upper, incr):
		sys.stderr.write(str(percent) + '\n')
		num_pick = int(round(percent/100. * num_events))
		for i in range(replicates):
			x = random.sample(event_results, num_pick)
			writer.writerow([percent, len(set(x)), num_pick])


