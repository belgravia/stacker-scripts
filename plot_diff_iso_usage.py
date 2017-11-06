import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.image as mpimg
import sys
import numpy as np

import matplotlib.patches as mplpatches 

plt.style.use('/pod/home/alison/mpl/BME163')

try:
    diffisofile = open(sys.argv[1])  # sirv consensus isoforms given uncorrected reads 
    plotname = sys.argv[2]
    if len(sys.argv) > 3:
        psl = open(sys.argv[3])
    else:
        psl = None
except:
    sys.stderr.write('usage: script.py diffisofile plotname [psl]\n')
    sys.exit(1)

def parse_psl(psl):
    info = []
    for line in psl:
        line = line.rstrip().split('\t')
        blocksizes = [int(n) for n in line[18].split(',')[:-1]]
        blockstarts = [int(n) for n in line[20].split(',')[:-1]]
        isoname = line[9][line[9].find('Isoform'):]
        if 'Isoform6' in isoname:
            continue
        info += [[blocksizes, blockstarts, isoname]] # [list of blocksizes, list of blockstarts, color]
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

def plot_blocks(data, panel, utr=False, height=.5, l=0.8):
    lower = sys.maxsize
    lim = 0
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
            isoname = line[2]
            panel.plot([min(starts)+2, max(starts)], [i]*2, 'k-', lw=l, color=colordict[isoname], alpha=0.9)
            for k in range(len(sizes)):  # each block of each read
                # if utr and line[2][k]:
                #     rectangle1 = mplpatches.Rectangle([starts[k], i-height/2], \
                #         sizes[k], height, facecolor='white', linewidth=0, zorder=11)
                #     rectangle2 = mplpatches.Rectangle([starts[k], i-.25/2],\
                #         sizes[k], 0.25, facecolor='black', linewidth=0, zorder=12)
                #     panel.add_patch(rectangle1)
                #     panel.add_patch(rectangle2)
                # else:
                rectangle = mplpatches.Rectangle([starts[k], i - height/2], \
                    sizes[k], height, facecolor=colordict[isoname], linewidth=0, alpha=0.9)
                panel.add_patch(rectangle)
                lim = max(lim, sizes[-1]+starts[-1])
                lower = min(lower, starts[0])
    panel.set_xlim(lower, lim + 100)
    panel.set_ylim(-1, len(data)+len(data)//10)

def plot_usage(panel, usage):
    counts, colors = [], []
    for i in usage.items():
        counts += [i[1]]  # count for this isoform
        colors += [colordict[i[0]]]  # color for this isoform
    theta = np.arange(0.0, 2*np.pi, 2*np.pi / len(usage))
    radii = [count / sum(counts) for count in counts]
    width = 2*np.pi/len(wt)
    bars = panel.bar(theta, radii, width=width, bottom=0.0, linewidth=0.5)
    for r, bar, col in zip(radii, bars, colors):
        bar.set_facecolor(col)
        bar.set_alpha(0.9)

colors_ = ['#4A3F3F',  '#FEB280', '#93B69C', '#5F94C6', '#F2D377']
colors_ = ['#8c510a', '#d8b365', '#FEB280', '#5F94C6', '#5ab4ac', '#01665e']
# retired colors '#A07068','#686644'
wt, mt, colordict = {}, {}, {}
for line in diffisofile:
    line = line.rstrip().split('\t')
    if 'Isoform6' in line[1]:
        continue
    wt[line[1]] = float(line[3])
    mt[line[1]] = float(line[4])
    colordict[line[1]] = colors_.pop()
# wt['Isoform4_E_S_2'] = 1
# mt['Isoform4_E_S_2'] = 7
print(wt)
print(mt)
smallestpie = 0.55
minsize = min(sum(list(wt.values())), sum(list(mt.values())))
wtsize = smallestpie * sum(list(wt.values())) / minsize
mtsize = min(smallestpie * sum(list(mt.values())) / minsize, 0.76)
print(wtsize, mtsize)
fig_1 = plt.figure(figsize=(10,7.5))
panel1 = plt.axes([-0.225 - (mtsize - smallestpie)/10, .02 + (mtsize - smallestpie)/2, 1, wtsize], frameon=False, polar=True)
panel2 = plt.axes([.225 + (mtsize - smallestpie)/10, .02 + (wtsize - smallestpie)/2, 1, mtsize], frameon=False, polar=True)
panel3 = plt.axes([0.01, 0.77, .99, 0.20], frameon=False)  # isoform blocks

plot_usage(panel1, wt)
plot_usage(panel2, mt)

panel1.set_xticklabels([])
panel1.set_yticklabels([])
panel1.grid(False)

panel2.set_xticklabels([])
panel2.set_yticklabels([])
panel2.grid(False)

panel3.set_xticklabels([])
panel3.set_yticklabels([])


plot_blocks(pack(parse_psl(psl)), panel3, True, l=1)

plt.savefig(plotname)
