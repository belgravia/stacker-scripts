import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import scipy.stats as stats
import matplotlib.image as mpimg
import sys

try:
    txt5 = open(sys.argv[1])
    txt3 = open(sys.argv[2])
    outfile = sys.argv[3]
except:
    sys.stderr.write('script.py sequences5.txt sequences3.txt outfile.png\n')
    sys.exit(1)

plt.style.use('/pod/home/alison/mpl/BME163')
path = '/pod/home/alison/mpl/'
A = mpimg.imread(path+'A.png')
C = mpimg.imread(path+'C.png')
T = mpimg.imread(path+'T.png')
G = mpimg.imread(path+'G.png')

fig_width, fig_height = 8, 3
panel_width, panel_height = 3, 1.5
x_lims = [-10, 10]
y_lims = [0, 2]
plt.figure(figsize=(fig_width,fig_height))
panel1 = plt.axes([0.15, 0.3, panel_width/fig_width, panel_height/fig_height])
panel2 = plt.axes([0.55, 0.3, panel_width/fig_width, panel_height/fig_height])

def revcomp(seq):
    seq = seq.replace('A', 'X').replace('C', 'Y')
    seq = seq.replace('T', 'A').replace('G', 'C')
    seq = seq.replace('X', 'T').replace('Y', 'G')
    return seq[::-1]

def update_counts(seq, strand):
    if strand in ['3', '-']:
        seq = revcomp(seq)
        for i in range(len(seq)):
            counts3[nuc_order.index(seq[i])][i] += 1
    else:
        for i in range(len(seq)):
            counts5[nuc_order.index(seq[i])][i] += 1
    if len(locations[chrom]) == 1:
        return sys.maxsize
    locations[chrom] = locations[chrom][1:]
    return locations[chrom][0][0]

def frequency(counts):
    temp = np.zeros(shape=[len(counts), len(counts[0])])
    for j in range(len(counts[0])):
        col_j = counts[:,j]
        temp[:,j] = col_j /float(sum(col_j))
    return temp

def stackheights(freq):
    heights = []
    for j in range(len(freq[0])):
        col = freq[:,j]
        heights += [2 - sum(col * np.log2(col)) * -1]
    return heights  # 1x20

def yheights(freq, stackheights):
    yh = np.zeros(shape=[4, 20])
    bh = np.zeros(shape=[4, 20])
    for j in range(len(freq[0])):
        col = freq[:,j]  # 4x1
        stackheight = stackheights[j]
        baseheights = stackheight * col
        bh[:,j] = baseheights
        for b1 in range(len(baseheights) - 1):
            for b2 in range(b1 + 1, len(baseheights)):  # pairwise comparisons
                if baseheights[b1] > baseheights[b2]:
                    yh[b1][j] += baseheights[b2]
                else:
                    yh[b2][j] += baseheights[b1]
    return yh, bh

nuc_order = ['A', 'C', 'T', 'G']
counts5 = np.zeros(shape=[4, 20], dtype=int)
counts3 = np.zeros(shape=[4, 20], dtype=int)

for line in txt5:
    line = line.rstrip()
    for i in range(len(line)):
        if line[i] != 'N':
            counts5[nuc_order.index(line[i])][i] += 1

for line in txt3:
    line = line.rstrip()
    for i in range(len(line)):
        if line[i] != 'N':
            counts3[nuc_order.index(line[i])][i] += 1

freq5 = frequency(counts5)
freq3 = frequency(counts3)

yh5, bh5 = yheights(freq5, stackheights(freq5))
yh3, bh3 = yheights(freq3, stackheights(freq3))
# print(bh3)

for i in np.arange(0, 20, 1):
    x = i - 10
    panel1.imshow(A, extent=[x, x+1, yh5[0][i], yh5[0][i] + bh5[0][i]], aspect='auto')
    panel1.imshow(T, extent=[x, x+1, yh5[2][i], yh5[2][i] + bh5[2][i]], aspect='auto')
    panel1.imshow(G, extent=[x, x+1, yh5[3][i], yh5[3][i] + bh5[3][i]], aspect='auto')
    panel1.imshow(C, extent=[x, x+1, yh5[1][i], yh5[1][i] + bh5[1][i]], aspect='auto')

    panel2.imshow(A, extent=[x, x+1, yh3[0][i], yh3[0][i] + bh3[0][i]], aspect='auto')
    panel2.imshow(T, extent=[x, x+1, yh3[2][i], yh3[2][i] + bh3[2][i]], aspect='auto')
    panel2.imshow(G, extent=[x, x+1, yh3[3][i], yh3[3][i] + bh3[3][i]], aspect='auto')
    panel2.imshow(C, extent=[x, x+1, yh3[1][i], yh3[1][i] + bh3[1][i]], aspect='auto')

panel1.tick_params(axis='both', which='both',\
                   bottom='on', labelbottom='on',\
                   left='on', labelleft='on',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')
panel2.tick_params(axis='both', which='both',\
                   bottom='on', labelbottom='on',\
                   left='off', labelleft='off',\
                   right='off', labelright='off',\
                   top='off', labeltop='off')

panel1.plot([0, 0], [0, 2], 'k-', lw=.5)
panel2.plot([0, 0], [0, 2], 'k-', lw=.5)

panel1.set_xlim(x_lims)
panel1.set_ylim(y_lims)
panel1.set_xlabel('Distance to\nSplice Site')
panel1.set_ylabel('Bits')
panel1.set_xticks(np.arange(-10, 15, 5))
panel1.set_title('5\'SS')

panel2.set_xlim(x_lims)
panel2.set_ylim(y_lims)
panel2.set_xlabel('Distance to\nSplice Site')
panel2.set_xticks(np.arange(-10, 15, 5))
panel2.set_title('3\'SS')

plt.savefig(outfile, dpi=1200)


# problems with code
# assumes that splice sites are 60 bp apart, at least
# if two bases have the same frequency, they will stack on top of each other
# revcomp is 7N time and i'm lazy