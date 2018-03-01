import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import sys
import scipy.stats as stats
from matplotlib import cm

# one off script

plt.style.use('/pod/home/alison/mpl/BME163')

try:
    txt = open(sys.argv[1])
    outfilename = sys.argv[2]
except:
    sys.stderr.write('usage: script.py filename out.base\n')
    sys.exit(1)

info = {}
for line in txt:
    line = line.rstrip().split('\t')
    if line[0] == '*':
        continue
    gene = line[0][:5]

    if gene not in info:
        info[gene] = {}
        info[gene]['names'] = []
        info[gene]['counts'] = []
    info[gene]['names'] += [line[0]]
    info[gene]['counts'] += [float(line[1])]


xaxis_pos = 1

for gene in info:
    plt.figure(figsize=(3,3))
    panel1 = plt.axes([0.15,0.15,2.25/3,2.25/3])

    counts = info[gene]['counts']
    cm_subsection = np.linspace(0.0, 1.0, len(counts)) 
    colors = [ cm.Set3(x) for x in cm_subsection ]
    plt.pie([c/float(sum(counts)) for c in counts], labels=info[gene]['names'], autopct='%1.1f%%', \
        shadow=False, colors=colors, pctdistance=0.8)

    panel1.tick_params(axis='both',which='both', \
                   bottom='off', labelbottom='off', \
                   left='off', labelleft='off', \
                   right='off', labelright='off', \
                   top='off', labeltop='off')


    plt.savefig(outfilename+'.'+gene+'.png')


# panel1.set_xlim(0,xlim)
# panel1.set_xlabel('cat')  # fm/ul
# panel1.set_xticks(np.arange(1, xlim))
# panel1.set_xticklabels(cats)

# ylim=5
# panel1.text(0.25, 0.85 * ylim, '$p$=' + pval, fontsize=6)
# panel1.set_ylim(-5,ylim)
# panel1.set_ylabel(r'log$_\mathregular{2}$(Fold change)')  # 

