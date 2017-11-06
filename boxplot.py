import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import sys
import scipy.stats as stats

# one off script

plt.style.use('/pod/home/alison/mpl/BME163')

try:
    txt = open(sys.argv[1])
    colnum0 = int(sys.argv[2])
    colnum1 = int(sys.argv[3])
    outfilename = sys.argv[4]
except:
    sys.stderr.write('usage: script.py filename colnum0 colnum1 out.png\n')
    sys.stderr.write('colnum0 contains categories like nmdpositive or nmdnegative, colnum1 contains values\n')
    sys.stderr.write('outfilename ending in png\n')
    sys.exit(1)

def spread(y, x_val, step=0.01, cutoff=0.05, bound=0.6):
    y.sort()
    placed_points = []
    for y_val in y:
        if len(placed_points) == 0:
            placed_points.append((x_val, y_val))
        else:
            potential_x_position = []
            for x_position in np.arange(x_val - bound, x_val + bound, step):
                distances = []
                for placed_point in placed_points:
                    distance = ((x_position - placed_point[0])**2 + \
                        (y_val*xlim/ylim - placed_point[1]*xlim/ylim)**2)**0.5
                    distances.append(distance)
                if min(distances) > cutoff:
                    potential_x_position.append(x_position)
            if len(potential_x_position) > 0:
                best_x_position=sorted(potential_x_position, key=lambda x: abs(x-x_val))[0]
                placed_points.append((best_x_position, y_val)) 
            else:
                # sys.stderr.write('point not placed: {}\n'.format(y_val))
                continue
    x = [pair[0] for pair in placed_points]
    y = [pair[1] for pair in placed_points]
    return x, y

values = {}
ylim = 0
ymin = 0
totalcounts = 0
for line in txt:  # read in psl info
    line = line.strip().split('\t')
    cat = line[colnum0]
    data = float(line[colnum1])
    totalcounts += data
    ylim = max(ylim, data)
    if not cat in values:
        values[cat] = []
    values[cat] += [data]
print('ylim', ylim)
ylim = 10
xlim = len(values.keys()) + 1
print('total counts: ' + str(totalcounts))
scalingfactor = (1/32+1/4+1+4)

plt.figure(figsize=(3,3))
panel1 = plt.axes([0.2,0.15,2.25/3,2.25/3])
xaxis_pos = 1
cats = sorted(values.keys())
for cat in cats:
    print(cat)
    normalized_exp = np.log2(values[cat])
    bp = panel1.boxplot(normalized_exp, \
        positions=[xaxis_pos], \
        patch_artist=True, widths=0.5)
    for box in bp['boxes']:
        box.set(edgecolor='black', facecolor='none', linewidth=0.75)
    for whisker in bp['whiskers']:
        whisker.set(color='black', linestyle='-', linewidth=0.75)
    for median in bp['medians']:
        median.set(color='black', linestyle='-', linewidth=0.75)
    for flier in bp['fliers']:
        flier.set(markersize=0)
    for cap in bp['caps']:
        cap.set(markersize=0, linewidth=.8)
    x, y = spread(normalized_exp, xaxis_pos, 0.03, 0.025, 0.42)
    panel1.scatter(x, y, s=2,\
                   facecolor=(0,0,0),\
                   edgecolor=(0,0,0),\
                   linewidth=0.0, alpha=0.3)
    # panel1.scatter(xaxis_pos, np.log2(float(cat)/scalingfactor*totalcounts+1), s=4,\
    #                 facecolor=(0,0,1),\
    #                 edgecolor=(0,0,1),\
    #                 linewidth=0.0)  # blue dots
    xaxis_pos += 1

panel1.tick_params(axis='both',which='both', \
                   bottom='on', labelbottom='on', \
                   left='on', labelleft='on', \
                   right='off', labelright='off', \
                   top='off', labeltop='off')

# print(stats.spearmanr([1/32] * len(values['0.03125']) + [1/4] * len(values['0.25']) +
#                         [1] * len(values['1']) + [4] * len(values['4']), 
#                         sorted(values['0.03125']) + sorted(values['0.25']) + 
#                         sorted(values['1']) + sorted(values['4'])))
# pearsonres = round(np.corrcoef([1/32] * len(values['0.03125']) + [1/4] * len(values['0.25']) +
#                         [1] * len(values['1']) + [4] * len(values['4']), 
#                         sorted(values['0.03125']) + sorted(values['0.25']) + 
#                         sorted(values['1']) + sorted(values['4']))[0][1], 3)
# panel1.text(1/4, 0.85 * ylim, '$r^2$=' + str(pearsonres), fontsize=8)


# wilcoxon = stats.wilcoxon(values[cats[0]], values[cats[1]])
# print(wilcoxon)


panel1.set_xlim(0,xlim)
panel1.set_xlabel('cat')  # fm/ul
panel1.set_xticks(np.arange(1, xlim))
panel1.set_xticklabels(cats)

ylim=5
panel1.text(0.25, 0.85 * ylim, '$p$=' + str('2.87e-5'), fontsize=6)
panel1.set_ylim(-5,ylim)
panel1.set_ylabel(r'log$_\mathregular{2}$(Fold change)')  # 

plt.savefig(outfilename)
