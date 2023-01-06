import os 
import sys 
import math
import pickle
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt

template = {'legend.numpoints': 1, 'axes.axisbelow': True, 'axes.labelcolor': '.15', 'ytick.major.size': 4.0, 'axes.grid': False, 'ytick.minor.size': 0.0, 'legend.scatterpoints': 1, 'axes.edgecolor': "black", 'grid.color': 'white', 'legend.frameon': False, 'ytick.color': '.15', 'xtick.major.size': 4.0, 'figure.facecolor': "#EAEAF2", 'xtick.color': '.15', 'xtick.minor.size': 3.0, 'xtick.direction': u'out', 'lines.solid_capstyle': u'round', 'grid.linestyle': u'-', 'image.cmap': u'Greys', 'axes.facecolor': "white", 'text.color': '.15', 'ytick.direction': u'out', 'axes.linewidth': 1.0,'xtick.major.width': 0.5, 'ytick.major.width': 0.5,}
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
    ax.tick_params(pad=2,length=2)
    sns.heatmap(matrix, cmap=cmap, linewidths=.1, linecolor="#AAAAAA", cbar=False, ax=ax) 
    ax.set_yticks(list(map(lambda x: x+0.5, range(16))))
    ax.set_yticklabels(targets,fontsize=12,rotation=0)
    ax.set_xticks(list(map(lambda x: x+0.5, range(78))))
    ax.set_xticklabels(orders,fontsize=12,rotation=90) 

def cbar(ax,cmap):
    bounds = np.linspace(0,1,6)
    cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, spacing='proportional', ticks=bounds)
    cb.ax.tick_params(labelsize=16,pad=2)
    #cb.set_ticks([0.5,2.5,4.5,,6.5,8.5,1.05])

def read_count_data(fname): 
    matrix  = []
    rawdata = [] 
    targets = [] 
    for line in open(fname):
        targets.append(line.rstrip().split(",")[0])
        row = list(map(float,line.rstrip().split(",")[1:]))
        matrix.append(list(map(lambda x: x*1.0/sum(row),row)))  
        rawdata.append(list(map(lambda x: x*1.0, row)))
    return matrix, rawdata, targets

if __name__ == "__main__":
    labels = []     
    with open(sys.argv[1]) as f:
        for line in f:
            labels.append(line.rstrip().split(",")[0]) 
                
    #Read count data
    matrix1,rd1, targets = read_count_data(sys.argv[2])
    others = set(labels) - set(targets)
    others = list(others)
    targets.sort()
    others.sort() 
    orders = targets + others
    
    fig  = plt.figure(figsize=(16,4)) 
    ax1  = fig.add_axes([0.1,0.1,0.8,0.8])
    ax2  = fig.add_axes([0.92,0.1,0.01,0.4]) 
    cmap = generate_cmap(["#000000","#227700","#55BB00","#77FF00"])
    heatmap(ax1,cmap,matrix1,0.0,1.0) 
    cbar(ax2,cmap)
    fig.savefig("{}.pdf".format(sys.argv[2].split("/")[-1].split(".")[0]),bbox_inches="tight")

    #Read count data
    matrix2, _, targets = read_count_data(sys.argv[3])
    fig  = plt.figure(figsize=(16,4)) 
    ax1  = fig.add_axes([0.1,0.1,0.8,0.8])
    ax2  = fig.add_axes([0.92,0.1,0.01,0.4]) 
    cmap = generate_cmap(["#000000","#227700","#55BB00","#77FF00"])
    heatmap(ax1,cmap,matrix2,0.0,1.0) 
    cbar(ax2,cmap)
    fig.savefig("{}.pdf".format(sys.argv[3].split("/")[-1].split(".")[0]),bbox_inches="tight")

    #Read count data
    matrix3, rd3, targets = read_count_data(sys.argv[4])
    fig  = plt.figure(figsize=(16,4)) 
    ax1  = fig.add_axes([0.1,0.1,0.8,0.8])
    ax2  = fig.add_axes([0.92,0.1,0.01,0.4]) 
    cmap = plt.cm.Blues
    heatmap(ax1,cmap,matrix3,0.0,1.0) 
    cbar(ax2,cmap)
    fig.savefig("{}.pdf".format(sys.argv[4].split("/")[-1].split(".")[0]),bbox_inches="tight")
    
    #Read count data
    matrix4, _, targets = read_count_data(sys.argv[5])
    fig  = plt.figure(figsize=(16,4)) 
    ax1  = fig.add_axes([0.1,0.1,0.8,0.8])
    ax2  = fig.add_axes([0.92,0.1,0.01,0.4]) 
    cmap = plt.cm.Blues
    heatmap(ax1,cmap,matrix4,0.0,1.0) 
    cbar(ax2,cmap)
    fig.savefig("{}.pdf".format(sys.argv[5].split("/")[-1].split(".")[0]),bbox_inches="tight")

    #Read count dat
    matrix5 = np.array(rd3) * 1.0/ (np.array(rd1)+1.0) 
    fig  = plt.figure(figsize=(16,4)) 
    ax1  = fig.add_axes([0.1,0.1,0.8,0.8])
    ax2  = fig.add_axes([0.92,0.1,0.01,0.4]) 
    cmap = plt.cm.Reds
    heatmap(ax1,cmap,matrix5,0.0,1.0) 
    cbar(ax2,cmap)
    fig.savefig("{}.pdf".format("atg_conversion"),bbox_inches="tight")

    #cmap        = plt.cm.Blues
    #cmaplist    = [cmap(i) for i in range(0,256)]
    #print len(cmaplist) 
    #cmaplist[0] = (1.0,1.0,1.0,1.0)
    #cmap        = cmap.from_list('Custom cmap', cmaplist, 100)
    
    #cmap        = plt.cm.Greens
    #cmaplist    = [cmap(i) for i in range(0,256,52)]
    #print len(cmaplist) 
    #cmaplist[0] = (0.0,0.0,0.0,1.0)
    #cmap        = cmap.from_list('Custom cmap', cmaplist)
