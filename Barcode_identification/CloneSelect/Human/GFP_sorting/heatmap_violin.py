import os 
import sys 
import math
import pickle
import random
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt

template = {'legend.numpoints': 1, 'axes.axisbelow': True, 'axes.labelcolor': '.15', 'ytick.major.size': 4.0, 'axes.grid': False, 'ytick.minor.size': 0.0, 'legend.scatterpoints': 1, 'axes.edgecolor': "black", 'grid.color': 'white', 'legend.frameon': False, 'ytick.color': '.15', 'xtick.major.size': 4.0, 'figure.facecolor': "#EAEAF2", 'xtick.color': '.15', 'xtick.minor.size': 3.0, 'xtick.direction': u'out', 'lines.solid_capstyle': u'round', 'grid.linestyle': u'-', 'image.cmap': u'Greys', 'axes.facecolor': "white", 'text.color': '.15', 'ytick.direction': u'out', 'axes.linewidth': 1.0,'xtick.major.width': 1.0, 'ytick.major.width': 1.0,}
sns.set(font = "Helvetica")
sns.set_context("poster", font_scale=1.2, rc={"lines.linewidth": 1.0}) 
sns.set_style(template)
sns.set_palette("colorblind")
colorblind = sns.color_palette("colorblind")

from matplotlib.colors import LinearSegmentedColormap
def generate_cmap(colors):
    color_list = []
    values = range(len(colors))
    vmax   = int(np.max(values)) 
    for v, c in enumerate(colors):
        color_list.append( (v*1.0/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)

def heatmap(ax,cmap,matrix,vmin,vmax): 
    ax = sns.heatmap(matrix, cmap=cmap, linewidths=.1, linecolor="#AAAAAA", cbar=False, ax=ax, vmin=0, vmax=0.8) 
    ax.tick_params(pad=2,length=3,labelsize=11)
    ax.set_yticks(list(map(lambda x: x+0.5, range(16))))
    ax.set_yticklabels(targets,rotation=0)
    ax.set_xticks([])

def cbar(ax,cmap):
    bounds = np.linspace(0,1.0,6)
    norm = matplotlib.colors.Normalize(0.0,1.0)  
    cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, spacing='proportional', ticks=bounds, orientation="horizontal", norm=norm)
    cb.ax.tick_params(labelsize=11,pad=2)
    cb.outline.set_linewidth(0.0)

def read_count_data(fname): 
    matrix  = []
    rawdata = [] 
    targets = [] 
    for line in open(fname):
        targets.append(line.rstrip().split(",")[0])
        row = list(map(float,line.rstrip().split(",")[1:]))
        matrix.append(list(map(lambda x: (x+1)*1.0/sum(row),row)))  
        rawdata.append(list(map(lambda x: x*1.0, row)))
    return matrix, rawdata, targets

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def custom_violin(ax, data, pos, widths=0.7, facecolors=[(1.0,1.0,1.0,1.0)] * 16, cc="#333333", bc="#111111"):
    parts = ax.violinplot(data, pos, vert=False, widths=widths, showmeans=False, showmedians=False, showextrema=False, bw_method=0.2)    
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(facecolors[i])
        pc.set_edgecolor((0,0,0,1))
        pc.set_alpha(1.0)
        pc.set_linewidth(0.8) 

if __name__ == "__main__":
    labels = []     
    with open(sys.argv[1]) as f:
        for line in f:
            labels.append(line.rstrip().split(",")[0])
    random.seed(0) 
    #Read count data
    matrix1,rd1, targets = read_count_data(sys.argv[2])
    others = set(labels) - set(targets)
    others = list(others)
    targets.sort()
    others.sort() 
    orders = targets + others
    
    fig  = plt.figure(figsize=(8,3)) 
    ax1  = fig.add_axes([0.1,0.1,0.8,0.8])
    #ax2  = fig.add_axes([0.91+0.2*0.5/9,0.06,0.2*8.0/9,0.03]) 
    ax2  = fig.add_axes([0.70,0.04,0.18,0.03]) 
    
    num = 0 
    for row in matrix1:
        positions = [random.random() * 0.5 for j in range(len(row))] 
        ax3 = fig.add_axes([0.93,0.85-num*0.05,0.18,0.05]) 
        ax3.zorder=num
        ax3.scatter(row,positions,s=4,facecolor=(0.0,0.0,0.0,0.5),linewidth=0.0,edgecolor=(0.5,0.5,0.5,1.0),zorder=100)
        combis = list(zip(row,positions))
        ax3.plot([combis[num][0],combis[num][0]],[-1,1],color="#0000FF",linewidth=2.0,alpha=1.0,zorder=100)
        num += 1
        for spine in ax3.spines.keys():
            spine = ax3.spines[spine]
            spine.set_visible(False) 
        ax3.set_xlim([0,1.0])
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax3.set_ylim([-0.1,0.6])
        if num % 2 == 0:
            ax3.patch.set_facecolor("#EEEEEE")
    
    ax3.tick_params(pad=2,length=4,labelsize=11)
    ax3.spines["bottom"].set_visible(True)
    ax3.spines["bottom"].set_position(("axes",-1.2))
    ax3.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
    cmap = generate_cmap(["#FFFFFF","#96FF96","#00FF00","#008000"])
    heatmap(ax1,cmap,matrix1,0.0,1.0) 
    cbar(ax2,cmap)
    fig.patch.set_alpha(0.0)
    fig.savefig("{}.pdf".format(sys.argv[2].split("/")[-1].split(".")[0]),bbox_inches="tight")

    
