import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.image as mpimg
import sys

import matplotlib.patches as mplpatches

fig_1 = plt.figure(figsize=(8, 5  ))#10,7.5))
panel1 = plt.axes([0.00, 0.79, 1, 0.25], frameon=True) 
panel2 = plt.axes([0.00, 0.01, 1, 0.77], frameon=True)

try:
    psl1 = open(sys.argv[1]) 
    gtf = open(sys.argv[2])  # SIRV_C_numbers.gtf
    gene = sys.argv[3]  # SIRV1, SIRV2, etc
    plotname = sys.argv[4] # sirv1.png
    if len(sys.argv) > 5:
        numreads = int(sys.argv[5])
    else:
        numreads = 50
except:
    sys.stderr.write('usage: python script.py psl1 gtf gene_to_plot plotname [numreads] \n')
    sys.exit(1)

plt.style.use('/pod/home/alison/mpl/BME163')

if gene == 'SIRV1':
    coord = (1000, 11700)
elif gene == 'SIRV2':
    coord = (1000, 6000)
elif gene == 'SIRV3':
    coord = (1000, 10000)
elif gene == 'SIRV4':
    coord = (1000, 15150)
elif gene == 'SIRV5':
    coord = (1000, 12000)
elif gene == 'SIRV6':
    coord = (1000, 12000)
elif gene == 'SIRV7':
    coord = (1000, 148000)
lim = coord[1] - coord[0]

def parse_psl(psl, color=False):
    info = []
    psl.readline()  # remove this line if there is no header in your psl
    for line in psl:
        line = line.rstrip().split('\t')
        if len(info) > numreads:
            break
        if gene in line[13]:
            blocksizes = [int(n) for n in line[18].split(',')[:-1]]
            blockstarts = [int(n) - coord[0] for n in line[20].split(',')[:-1]]
            if color:
                if 'p_p' in line[9]:
                    col = (0,0,0)
                else:
                    continue
                # elif 'p_n' in line[9]:  # uncomment these
                #     col = (0,1,0)
                # elif 'n_p' in line[9]:
                #     col = (0,0,1)
                # elif 'n_n' in line[9]:
                #     col = (1,0,0)
                # else:
                #     sys.stderr.write('Perhaps do not call this with the color option')
                #     continue
                info += [[blocksizes, blockstarts, col]]
            else:
                info += [[blocksizes, blockstarts]] # [list of blocksizes, list of blockstarts, color]
    return info

def parse_gtf(gtf):
    info = []
    current_transcript = ''
    for line in gtf:
        line = line.rstrip().split('\t')
        if len(line) < 8:
            continue
        start = line[8].find('transcript_id')
        if start < 0:  # not a transcript
            continue
        transcript_id = line[8][start+len('transcript_id')+2:]
        transcript_id = transcript_id[:transcript_id.find('";')]
        if gene in transcript_id:
            if transcript_id != current_transcript:
                current_transcript = transcript_id
                info += [[[],[],[]]]  # size_list, start_list, utr_boolean
            if line[2] == 'exon':
                info[-1][2] += [0]  # boolean
            elif line[2] == 'UTR':
                info[-1][2] += [1]  # will be plotted thinner later
            else:
                continue  # didn't take any other types into consideration
            info[-1][0] += [int(line[4]) - int(line[3])]
            info[-1][1] += [int(line[3]) - coord[0]]
    return info

def pack(data, rev=True, color=False):
    starts = [max(d[1]) for d in data] if rev else [min(d[1]) for d in data] # sort by right or left end
    data = [d for (s,d) in sorted(zip(starts, data))]
    packed = [[ data[0] ]]
    ends = [max(data[0][1]) + data[0][0][data[0][1].index(max(data[0][1]))]]
    for i in range(1, len(data)):
        min_start = min(data[i][1])
        end = max(data[i][1]) + data[i][0][data[i][1].index(max(data[i][1]))]
        pos = -1
        for j in range(len(packed)):
            if ends[j] < min_start:
                pos = j  # pack this read with the read at position j
                break
        if pos >= 0:
            packed[pos] += [data[i]]
            ends[pos] = end
        else:
            packed += [[data[i]]]
            ends += [end]
    return packed

def plot_blocks(data, panel, utr=False, height=.5, l=0.8, color=False):
    panel.set_xlim(1, lim)
    panel.set_ylim(-1, len(data)+1)
    panel.tick_params(axis='both', which='both',\
                       bottom='off', labelbottom='off',\
                       left='off', labelleft='off',\
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    for i in range(len(data)):  # each line
        read = data[i]
        for j in range(len(read)):  # each read on each line
            line = read[j]
            sizes = line[0]
            starts = line[1]
            panel.plot([min(starts)+2, max(starts)], [i]*2, 'k-', lw=l)
            for k in range(len(sizes)):  # each block of each read
                if utr and line[2][k]:
                    rectangle1 = mplpatches.Rectangle([starts[k], i-height/2], \
                        sizes[k], height, facecolor='white', linewidth=0, zorder=11)
                    rectangle2 = mplpatches.Rectangle([starts[k], i-.25/2],\
                        sizes[k], 0.25, facecolor='black', linewidth=0, zorder=12)
                    panel.add_patch(rectangle1)
                    panel.add_patch(rectangle2)
                elif color:
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2], \
                        sizes[k], height, facecolor=line[2], linewidth=0)
                    panel.add_patch(rectangle)                   
                else:
                    rectangle = mplpatches.Rectangle([starts[k], i - height/2], \
                        sizes[k], height, facecolor='black', linewidth=0)
                    panel.add_patch(rectangle)    


plot_blocks(pack(parse_gtf(gtf)), panel1, True, l=1)
plot_blocks(pack(parse_psl(psl1, color=False)), panel2, l=0.08, color=False)

plt.savefig(plotname)
