import sys, csv
from subprocess import call

try:
    n = int(sys.argv[1])
    command = sys.argv[2].split(',')
    numeach = int((len(sys.argv) - 3)/n)
except:
    sys.stderr.write('usage: thisscript.py n script_to_run.py args0* args1* ... argsn*\n')
    sys.stderr.write('script_to_run.py args00 args10 argsn0\n')
    sys.exit(1)

argsets = []
for i in range(n):
	argsets += [sys.argv[3+i*numeach:3+(i+1)*numeach]]

for i in range(numeach):
	args = []
	for j in range(n):
		args += [argsets[j][i]]
	# subprocess.call(command + args + ['> {}.stdout 2> {}.stderr'.format(i, i)])
	print(args[0])
	call(['python', '/pod/pstore/groups/brookslab/bin/NanoSim_Wrapper.py', '--count', args[0], \
		'--fasta_file', '~/bl/atang/reference/gencode.v24.transcripts.fa', \
		'--model_prefix', '~/bl/kevynhart/NanoSim/Hg38_Model/Hg38_training', \
		'--outdir', 'ns.sf3b1.hg38.'+str(i),'--combined', 'simulated.sf3b1.fasta'+str(i)])