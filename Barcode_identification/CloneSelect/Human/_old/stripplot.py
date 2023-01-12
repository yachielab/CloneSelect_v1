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

random.seed(0) 
template = {'legend.numpoints': 1, 'axes.axisbelow': True, 'axes.labelcolor': '.15', 'ytick.major.size': 4.0, 'axes.grid': False, 'ytick.minor.size': 0.0, 'legend.scatterpoints': 1, 'axes.edgecolor': "black", 'grid.color': 'white', 'legend.frameon': False, 'ytick.color': '.15', 'xtick.major.size': 4.0, 'figure.facecolor': "#EAEAF2", 'xtick.color': '.15', 'xtick.minor.size': 3.0, 'xtick.direction': u'out', 'lines.solid_capstyle': u'round', 'grid.linestyle': u'-', 'image.cmap': u'Greys', 'axes.facecolor': "white", 'text.color': '.15', 'ytick.direction': u'out', 'axes.linewidth': 1.0,'xtick.major.width': 1.0, 'ytick.major.width': 1.0,}
sns.set(font = "Helvetica")
sns.set_context("poster", font_scale=1.2, rc={"lines.linewidth": 1.0}) 
sns.set_style(template)
sns.set_palette("colorblind")
colorblind = sns.color_palette("colorblind")

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

    fig = plt.figure(figsize=(16,4)) 
    ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
    ax1.tick_params(length=6.0, pad=4, labelsize=16)
    ax1.tick_params(axis="x", length=0.0, pad=4, labelsize=16)
    num = 0 
    for row in matrix1:
        positions = [((random.random() * 0.70) - 0.4) + num for j in range(len(row))] 
        ax1.scatter(positions,row,s=40,c="#AAAAAA",linewidth=0)
        combis = list(zip(row,positions))
        #combis.sort() 
        ax1.scatter([combis[num][1]],[combis[num][0]],s=40,c=colorblind[0])
        num += 1
    ax1.spines["top"].set_visible(False) 
    ax1.spines["bottom"].set_visible(False) 
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_bounds(0, 1)
    ax1.set_xlim(left=-0.8,right=len(targets)-0.2)
    ax1.set_xticks(list(range(len(targets)))) 
    ax1.set_xticklabels(targets,rotation=90)
    ax1.set_ylim(-0.05,1.05) 
    ax1.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    fig.patch.set_alpha(0.0)
    fig.savefig("{}_scatter.pdf".format(sys.argv[2].split("/")[-1].split(".")[0]),bbox_inches="tight")


