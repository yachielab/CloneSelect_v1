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
    matrix1, rd1, targets = read_count_data(sys.argv[2])
    matrix2, rd2, targets = read_count_data(sys.argv[3])

    n = 0 
    for row1, row2 in zip(matrix1,matrix2):
        row1.remove(row1[n]) 
        row2.remove(row2[n])

