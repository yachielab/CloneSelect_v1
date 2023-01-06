import os 
import re
import sys
import math
import logomaker
import Levenshtein
import pandas as pd
import itertools
import Levenshtein
import collections
import multiprocessing as mp 
import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm
from Bio import pairwise2
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ["Helvetica","Arial","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['figure.figsize']  = [4,4]
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams["axes.labelcolor"] = "#000000"
matplotlib.rcParams["axes.linewidth"]  = 1.0 
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0
cmap1 = plt.cm.tab10
cmap2 = plt.cm.Set3  
colors1 = [cmap1(i) for i in range(0,10)]
colors2 = [cmap2(i) for i in range(0,12)] 
def calc_bitscore(fname, start, end, strand): 
    pos_count_dict = {}
    record = SeqIO.read(fname,'abi')
    j = 0 
    
    if strand == "Forward":
        for i, pos in enumerate(record.annotations['abif_raw']["PLOC1"]):
            if start-1 <= i < end:
                pos_count_dict[j] = {} 
                pos_count_dict[j]["G"] = record.annotations['abif_raw']["DATA9"][pos] 
                pos_count_dict[j]["A"] = record.annotations['abif_raw']["DATA10"][pos]
                pos_count_dict[j]["T"] = record.annotations['abif_raw']["DATA11"][pos]  
                pos_count_dict[j]["C"] = record.annotations['abif_raw']["DATA12"][pos] 
                j += 1
    
    elif strand == "Reverse":
        for i, pos in enumerate(list(reversed(record.annotations['abif_raw']["PLOC1"]))):
            if start-1 <= i < end:
                pos_count_dict[j] = {} 
                pos_count_dict[j]["C"] = record.annotations['abif_raw']["DATA9"][pos] 
                pos_count_dict[j]["T"] = record.annotations['abif_raw']["DATA10"][pos]
                pos_count_dict[j]["A"] = record.annotations['abif_raw']["DATA11"][pos]  
                pos_count_dict[j]["G"] = record.annotations['abif_raw']["DATA12"][pos] 
                j += 1

    for pos in pos_count_dict:
        v = 100000 / sum(list(pos_count_dict[pos].values()))
        pos_count_dict[pos]["G"] = int(v*pos_count_dict[pos]["G"])
        pos_count_dict[pos]["A"] = int(v*pos_count_dict[pos]["A"])
        pos_count_dict[pos]["T"] = int(v*pos_count_dict[pos]["T"])
        pos_count_dict[pos]["C"] = int(v*pos_count_dict[pos]["C"])
        
        while sum(list(pos_count_dict[pos].values())) < 100000:
            for key in pos_count_dict[pos]:
                pos_count_dict[pos][key] += 1
    
    outname = "count_matrix"
    with open(outname + ".txt","w") as o:
        positions = list(pos_count_dict.keys()) 
        positions.sort()
        o.write("A\tC\tG\tT\n") 
        for pos in positions:
            values = [] 
            for nucl in ["A", "C", "G", "T"]:
                values.append(pos_count_dict[pos][nucl]+1) 
            o.write("\t".join(list(map(str,values))) + "\n")  
    
    df  = pd.read_csv(outname + ".txt",sep="\t")
    return df

def sanger_to_logo(fname, start, end, strand, ax, outname=None):
    pos_count_dict = {}
    record = SeqIO.read(fname,'abi')
    j = 0 
    
    if strand == "Forward":
        for i, pos in enumerate(record.annotations['abif_raw']["PLOC1"]):
            if start-1 <= i < end:
                pos_count_dict[j] = {} 
                pos_count_dict[j]["G"] = record.annotations['abif_raw']["DATA9"][pos] 
                pos_count_dict[j]["A"] = record.annotations['abif_raw']["DATA10"][pos]
                pos_count_dict[j]["T"] = record.annotations['abif_raw']["DATA11"][pos]  
                pos_count_dict[j]["C"] = record.annotations['abif_raw']["DATA12"][pos] 
                j += 1
    
    elif strand == "Reverse":
        for i, pos in enumerate(list(reversed(record.annotations['abif_raw']["PLOC1"]))):
            if start-1 <= i < end:
                pos_count_dict[j] = {} 
                pos_count_dict[j]["C"] = record.annotations['abif_raw']["DATA9"][pos] 
                pos_count_dict[j]["T"] = record.annotations['abif_raw']["DATA10"][pos]
                pos_count_dict[j]["A"] = record.annotations['abif_raw']["DATA11"][pos]  
                pos_count_dict[j]["G"] = record.annotations['abif_raw']["DATA12"][pos] 
                j += 1

    for pos in pos_count_dict:
        v = 100000 / sum(list(pos_count_dict[pos].values()))
        pos_count_dict[pos]["G"] = int(v*pos_count_dict[pos]["G"])
        pos_count_dict[pos]["A"] = int(v*pos_count_dict[pos]["A"])
        pos_count_dict[pos]["T"] = int(v*pos_count_dict[pos]["T"])
        pos_count_dict[pos]["C"] = int(v*pos_count_dict[pos]["C"])
        
        while sum(list(pos_count_dict[pos].values())) < 100000:
            for key in pos_count_dict[pos]:
                pos_count_dict[pos][key] += 1
    
    if outname == None:
        outname = "count_matrix"
    else:
        pass
    
    with open(outname + ".txt","w") as o:
        positions = list(pos_count_dict.keys()) 
        positions.sort()
        o.write("A\tC\tG\tT\n") 
        for pos in positions:
            values = [] 
            for nucl in ["A", "C", "G", "T"]:
                values.append(pos_count_dict[pos][nucl]) 
            o.write("\t".join(list(map(str,values))) + "\n")  
    
    df  = pd.read_csv(outname + ".txt",sep="\t")

    logo = logomaker.Logo(df,
            font_name='Helvetica',
            color_scheme='classic',
            vpad=.0,
            width=.8,
            ax=ax)
    return ax, df

if __name__ == "__main__":
    hamming = lambda x,y : sum([1 for a,b in zip(x,y) if a == b or (a == "G" and b == "A")])
    intensity_dict = {} 
    with open(sys.argv[1]) as f:
        f.readline() 
        for line in f:
            line = line.rstrip().split("\t")
            intensity_dict[line[0]] = {}
            intensity_dict[line[0]]["Sanger"]  = {}
            intensity_dict[line[0]]["Sanger"]["file_name"] = line[1]
            intensity_dict[line[0]]["Barcode"] = line[2]
            intensity_dict[line[0]]["Values"]  = {}
            intensity_dict[line[0]]["Values"]["Rep1"] = line[3]
            intensity_dict[line[0]]["Values"]["Rep2"] = line[4]
            intensity_dict[line[0]]["Values"]["Rep3"] = line[5]

    common  = "TCAACTGCTAGCATCCGCCG"
    common2 = "ATGGC" 
    #for key in ["91_pool_100_BC21#4"]:
    for key in intensity_dict:
        if intensity_dict[key]["Sanger"]["file_name"] != "(N.A.)":
            fname  = "./pool_1580/"+ intensity_dict[key]["Sanger"]["file_name"]
            record = SeqIO.read(fname,'abi')
            seq = ""
            for i, pos in enumerate(record.annotations['abif_raw']["PLOC1"]):
                nucls  = ["G","A","T","C"] 
                values = [record.annotations['abif_raw']["DATA{}".format(j)][pos] for j in range(9,13)]
                seq += nucls[values.index(max(values))]
            
            
            values  = [] 
            seq     = seq.translate(str.maketrans("ATGCatgc","TACGtacg"))[::-1]
            barcode = intensity_dict[key]["Barcode"]
            for i in range(len(seq) - 20):
                values.append(hamming(seq[i:i+20],common))
            
            index = values.index(max(values)) + 20
            start = index + 1
            if common2 in seq:
                end = seq.find(common2) + 3
            else:
                end = index + 20
            alignments = pairwise2.align.globalms(barcode,seq[start-1:end], 2, 0, -1, -1)
            query   = alignments[0][0] 
            subject = alignments[0][1]
            #print(alignments[0])
            df   = calc_bitscore(fname, start, end, "Reverse")
            bits = []
            for c, index in zip(query,df.index.values):
                #print(df.loc[index].sum())
                bit = math.log2(4) + sum(list(map(lambda x: (x/df.loc[index].sum()) * math.log2(x/df.loc[index].sum()), df.loc[index])))
                if c == "G":
                    pass
                else:
                    bits.append(bit) 
                #print(bit * 1/0.69, df.loc[index])
            mean = np.mean(bits)
            ratio = df.loc[len(query.replace("-",""))-3]["A"]/df.loc[len(query.replace("-",""))-3].sum()
            print(key,query,subject,*intensity_dict[key]["Values"].values(),1-hamming(query,subject)/len(query),mean,ratio,sep="\t")
            """
            fig   = plt.figure(figsize=(0.25,1))
            ax    = fig.add_axes([0.1, 0.1, 0.8*(end-start+1), 0.75])
            ax, df = sanger_to_logo(fname, start, end, "Reverse", ax, outname=None)
            ax.set_yticks([])
            ax.spines["top"].set_visible(False) 
            ax.spines["right"].set_visible(False) 
            ax.spines["left"].set_visible(False) 
            fig.savefig("logo/{}".format(key + ".pdf"), bbox_inches="tight")
            """
            plt.close()
        
        else:
            print(key,"N.A.","N.A.",*intensity_dict[key]["Values"].values(),"N.A.","N.A.","N.A.",sep="\t")

