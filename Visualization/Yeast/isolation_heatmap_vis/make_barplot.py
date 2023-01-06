import os 
import re
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import patchworklib as pw
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm
from matplotlib import font_manager
import matplotlib.ticker as ticker
import seaborn as sns 

#font_manager.fontManager.addfont("/Library/Fonts/Helvetica.ttf")
matplotlib.rcParams['ps.fonttype']       = 42
matplotlib.rcParams['pdf.fonttype']      = 42
matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.sans-serif']   = ["Arial", "Helvetica","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['figure.figsize']    = [4,4]
matplotlib.rcParams['font.size']         = 20
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.25
matplotlib.rcParams["xtick.major.width"] = 1.25
matplotlib.rcParams["xtick.minor.width"] = 1.25
matplotlib.rcParams["ytick.major.width"] = 1.25
matplotlib.rcParams["ytick.minor.width"] = 1.25
matplotlib.rcParams['xtick.major.pad']   = 4
matplotlib.rcParams['ytick.major.pad']   = 6
matplotlib.rcParams['xtick.major.size']  = 6
matplotlib.rcParams['ytick.major.size']  = 6
matplotlib.rcParams['xtick.minor.size']  = 4
matplotlib.rcParams['ytick.minor.size']  = 4

def generate_cmap(colors):
    color_list = []
    values = range(len(colors))
    vmax   = int(np.max(values)) 
    for v, c in enumerate(colors):
        color_list.append( (v*1.0/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)

cmap1 = plt.cm.tab10
cmap2 = plt.cm.Set3  
colors1 = [cmap1(i) for i in range(0,10)]
colors2 = [cmap2(i) for i in range(0,12)] 
colorblind = ["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"]
muted      = ["#4878D0", "#EE854A", "#6ACC64", "#D65F5F", "#956CB4", "#8C613C", "#DC7EC0", "#797979", "#D5BB67", "#82C6E2"]
bright     = ["#023EFF", "#FF7C00", "#1AC938", "#E8000B", "#8B2BE2", "#9F4800", "#F14CC1", "#A3A3A3", "#FFC400", "#00D7FF"]
pastel     = ["#A1C9F4", "#FFB482", "#8DE5A1", "#FF9F9B", "#D0BBFF", "#DEBB9B", "#FAB0E4", "#CFCFCF", "#FFFEA3", "#B9F2F0"]
okabe_ito  = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

Set1       = sns.color_palette('Set1'   , 8).as_hex()
Set2       = sns.color_palette('Set2'   , 8).as_hex()
Set3       = sns.color_palette('Set3'   , 8).as_hex()
accent     = sns.color_palette('Accent' , 8).as_hex()
Dark2      = sns.color_palette('Dark2'  , 8).as_hex()

ax  = pw.Brick(figsize=(1, 4.0))
ax1 = pw.Brick(figsize=(1, 0.4))
ax2 = pw.Brick(figsize=(1, 0.4))
ax3 = pw.Brick(figsize=(1, 0.4))
ax4 = pw.Brick(figsize=(1, 0.4))

cmap1 = plt.cm.YlGn
cmap1 = [cmap1(i) for i in range(0,256)][0:]
cmap2 = plt.cm.PuBu
cmap2 = [cmap2(i) for i in range(0,256)][0:]
cmap3 = plt.cm.Wistia
cmap3 = [cmap3(i) for i in range(0,256)][0:]

#cmap4 = plt.cm.gist_rainbow
if sys.argv[1] == "100":
    max_lcs = 26
    cmap4 = generate_cmap(Set2) #For pool100 

elif sys.argv[1] == "1580":
    max_lcs = 33 
    cmap4   = generate_cmap(pastel) #For pool1580

cmap4 = [cmap4(i) for i in range(0,256)][0:]

if __name__ == "__main__":
    targetseq_cloneid_dict = {} 
    with open(sys.argv[2]) as f: #load cloneid_targetseq.tsv
        f.readline() 
        for line in f:
            line = line.rstrip().split("\t") 
            targetseq_cloneid_dict[line[1].replace("-","")] = line[0] 

    with open(sys.argv[3]) as f:
        f.readline() 
        labels = [] 
        lcs    = [] 
        pos    = 0 
        ints   = []
        bits   = []
        lines  = f.readlines() 
        
        _ids      = []  
        new_lines = [] 
        for line in lines:
            line       = line.rstrip().split("\t")
            target_seq = line[1].replace("-","")
            try:
                line.append(targetseq_cloneid_dict[target_seq]) 
                _ids.append(int(line[-1].split(" ")[1]))   
            except:
                line.append("Clone 100000000") 
            new_lines.append(line) 
            
        for line in new_lines:
            lc = int(line[0].split("-")[0][2:]) if line[0][0:2] == "BC" else -1
            if line[0][0:2] == "BC":
                lcs.append(int(lc))
            else:
                lcs.append(0)

        lcs_count_dict = {} 
        for lc in set(lcs):
            lcs_count_dict[lc] = lcs.count(lc) 

        count         = 1
        pre_count     = 1 
        pre_pre_count = 1 
        pre_d         = -1
        new_lines.sort(key=lambda x: int(x[-1].split(" ")[-1])) 
        
        for line in new_lines: 
            _id = line[-1]
            labels.append(line[0])
            lc = int(line[0].split("-")[0][2:]) if line[0][0:2] == "BC" else 0
            values = list(map(float,line[3:6]))
            if line[1] == "N.A.":
                if line[0] == "RFP-negative":
                    ax.bar([pos],[np.mean(values)],width=0.7,facecolor="#BAB0AC",zorder=50)
                elif line[0] == "RFP-positive":
                    ax.bar([pos],[np.mean(values)],width=0.7,facecolor="#FF515E",zorder=50)
                elif line[0] == "Blank, medium only [SC-HLU+Ade]":
                    ax.bar([pos],[np.mean(values)],width=0.7,facecolor="#FFE699",zorder=50)

                ax.scatter([pos-0.2,pos,pos+0.2], values, s=10, lw=0.8, facecolor="#30303000", edgecolor="k", zorder=100)  
                ax1.bar([pos], [1], width=1.0, facecolor = "#FFFFFF", lw=0.5, edgecolor="#AAAAAA")
                ax1.plot([pos-0.5,pos+0.5],[1,0],lw=1,c="k") 
                ax2.bar([pos], [1], width=1.0, facecolor = "#FFFFFF", lw=0.5, edgecolor="#AAAAAA")
                ax2.plot([pos-0.5,pos+0.5],[1,0],lw=1,c="k") 
                ax3.bar([pos], [1], width=1.0, facecolor = "#FFFFFF", lw=0.5, edgecolor="#AAAAAA")
                ax3.plot([pos-0.5,pos+0.5],[1,0],lw=1,c="k") 

            else:
                ax.bar([pos],[np.mean(values)],width=0.7,facecolor="#FF88BD")
                ax.scatter([pos-0.2, pos, pos+0.2], values, s=10, lw=0.8, facecolor="#30303000", edgecolor="k", zorder=100) 
                a = int(255 * (float(line[6]) / 0.5))
                b = int(255 * ((float(line[7])-0.8) / 1.2))
                c = int(255 * ((float(line[8])-0.5) / 0.5))
                ax1.bar([pos], [1], width=1.0, facecolor = cmap1[a if a < 256 else 255], lw=0.5, edgecolor="#CCCCCC")
                ax2.bar([pos], [1], width=1.0, facecolor = cmap2[b if b < 256 else 255], lw=0.5, edgecolor="#CCCCCC")
                ax3.bar([pos], [1], width=1.0, facecolor = cmap3[c if c < 256 else 255], lw=0.5, edgecolor="#cccccc")
                ints.append(np.mean(values)) 
                bits.append(float(line[7]))
          
            d = int(255 * (lc / max_lcs)) 
            if lc > 0 and d != pre_d:
                ax4.bar([pos-0.5], [1], width=lcs_count_dict[lc], facecolor=cmap4[d if d < 256 else 255], lw=0.5, edgecolor="k", align="edge")
                if lcs_count_dict[lc] == 1:
                    if count == 1:
                        ax4.text(((pos-0.5+pos-0.5+lcs_count_dict[lc])/2)-0.01/len(lines)+0.7, -0.3, "{}".format(_id), ha="right", va="top", color="k", rotation=60) 
                    else:
                        ax4.text(((pos-0.5+pos-0.5+lcs_count_dict[lc])/2)-0.01/len(lines)+0.7-0.2, -0.3, "{}".format(_id), ha="right", va="top", color="k", rotation=60) 
                else:
                    ax4.text(((pos-0.5+pos-0.5+lcs_count_dict[lc])/2)-0.01/len(lines) + 0.5, -0.3, "{}".format(_id), ha="right", va="top", color="k", rotation=60) 
                pre_pre_count = pre_count 
                pre_count = count 
                count = 1
            else:
                count += 1
                pass
            pre_d = d
            pos += 1    

        ax.set_xticks([])
        ax.set_yscale("log") 
        ax.set_ylim(100,200000) 
        ax.spines["top"].set_visible(False) 
        ax.spines["right"].set_visible(False) 
        ax1.set_xticks([]) 
        ax1.set_yticks([]) 
        ax2.set_xticks([]) 
        ax2.set_yticks([]) 
        ax3.set_xticks([]) 
        ax3.set_yticks([]) 
        ax4.set_xticks([]) 
        ax4.set_yticks([]) 

        ax.set_xlim(-0.5,pos-0.5) 
        ax1.set_xlim(-0.5,pos-0.5) 
        ax2.set_xlim(-0.5,pos-0.5) 
        ax3.set_xlim(-0.5,pos-0.5) 
        ax4.set_xlim(-0.5,pos-0.5) 
        ax1.set_ylim(0,1.0)
        ax2.set_ylim(0,1.0)
        ax3.set_ylim(0,1.0)
        ax4.set_ylim(0,1.0) 
        pw.Brick._figure.patch.set_alpha(0.0)
        
        pw.param["margin"] = None
        ax.change_aspectratio((pos*4.2/20,  2.5)) 
        ax1.change_aspectratio((pos*4.2/20, 0.25)) 
        ax2.change_aspectratio((pos*4.2/20, 0.25)) 
        ax3.change_aspectratio((pos*4.2/20, 0.25)) 
        ax4.change_aspectratio((pos*4.2/20, 0.175)) 
        
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0,numticks=100))
        ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0,subs=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), numticks=100)) 
        
        fig = (ax/ax1/ax2/ax3/ax4)
        fig.set_supylabel("Norm. RFP intensity (a.u.)") 
       
        print(sys.argv[3]) 
        fig.savefig("{}.pdf".format(sys.argv[3].split("/")[-1].split(".")[0]), bbox_inches="tight") 
        print("finish") 
