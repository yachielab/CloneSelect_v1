import os 
import re
import sys
import time 
import argparse
import Levenshtein
import collections
import subprocess 
import multiprocessing as mp 
import pickle
import math
import numpy as np
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ["Helvetica","Arial","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['figure.figsize']  = [4,3]
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams["axes.labelcolor"] = "#000000"
matplotlib.rcParams["axes.linewidth"]  = 1.0 
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0
cmap1 = plt.cm.tab10
cmap2 = plt.cm.Set3  
colors1 = [cmap1(i) for i in range(0,10)]
colors2 = [cmap2(i) for i in range(0,12)] 

def read_fastq(fastq_name):
    seq_dict = {} 
    with open(fastq_name) as f:
        n = 0 
        for line in f:
            if line[0] == "@" and n % 4 == 0:
                key = line[1:].rstrip() 
                key = key.split(" ")[0] 
                key = key.replace(":","_")
                seq_dict[key] = [line[1:].rstrip()] 
            else:
                seq_dict[key].append(line.rstrip()) 
            n += 1
    return seq_dict 

class SAMPLE(object):
    def __init__(self):
        pass
    
    def __getitem__(self, key):
        return self.__dict__[key]
    
    def __setitem__(self, key, value):
        self.__dict__[key] = value 

    def merge(self):   
        R1_dict = read_fastq("R1.fastq") 
        R2_dict = read_fastq("R2.fastq") 
        with open("R1_tmp.fastq","w") as R1_tmp, open("R2_tmp.fastq","w") as R2_tmp:
            for key in list(set(R1_dict.keys()) & set(R2_dict.keys())):
                for i in range(4):
                    R1_tmp.write(R1_dict[key][i] + "\n") 
                    R2_tmp.write(R2_dict[key][i] + "\n") 
        os.system("flash --max-overlap 300 --allow-outies --max-mismatch-density 0.3 {} {} > merge_result.txt".format("R1_tmp.fastq","R2_tmp.fastq"))
        os.system("mv out.extendedFrags.fastq merged.fastq")
        os.system("cat merged.fastq | awk \'NR % 4 == 1 {{split($0, a, \" \"); gsub(\":\",\"_\",a[1]); print \">\" substr(a[1],2)}} NR % 4 == 2 {{print $0}}\' > merged.fasta")
        os.system("mv out.notCombined_1.fastq notMerged_1.fastq")
        os.system("mv out.notCombined_2.fastq notMerged_2.fastq")
        os.system("rm R1_tmp.fastq R2_tmp.fastq") 
     
    def alignment(self, target, method="Needle", num_threads=1): 
        if "reference" not in os.listdir("."):
            os.mkdir("reference")
        else:
            pass
        
        with open("reference/reference.fasta","w") as o:
            o.write(">reference\n")
            o.write(self.references[0] + "\n") 
         
        if method == "Needle" and num_threads == 1:
            if target == "R2":
                os.system("needleall ./reference/reference.fasta {}.fasta -sreverse2 True -gapopen 20.0 -gapextend 2.0 -endopen 0.0 -endextend 0.0 -aformat markx3 -outfile {}.needle".format(target,target)) 
            else:
                os.system("needleall ./reference/reference.fasta {}.fasta -gapopen 20.0 -gapextend 2.0 -endopen 0.0 -endextend 0.0 -aformat markx3 -outfile {}.needle".format(target,target)) 
        
        if method == "Needle" and num_threads > 1:
            commands  = []
            num_lines = sum([1 for line in open("{}.fasta".format(target))])
            os.system("split -l {} {}.fasta tmp_".format(2 * (int((num_lines/2)/num_threads) + 1), target))
            for fasta in [f for f in os.listdir(".") if "tmp_" in f]:
                if target == "R2":
                    commands.append("needleall ./reference/reference.fasta {} -sreverse2 True -gapopen 20.0 -gapextend 2.0 -endopen 0.0 -endextend 0.0 -aformat markx3 -outfile {}.needle".format(fasta,fasta)) 
                else:
                    commands.append("needleall ./reference/reference.fasta {} -gapopen 20.0 -gapextend 2.0 -endopen 0.0 -endextend 0.0 -aformat markx3 -outfile {}.needle".format(fasta,fasta))
            
            os.system("&".join(commands))
            os.system("cat tmp_*.needle > {}.needle".format(target))
            os.system("rm tmp*")    
          
    def extract_target_regions(self, target, method="Needle", fname="target_region", lim_quality=-50):
        left_margin  = 0
        right_margin = len(self.references[0])
        min_match    = self.min_match
        
        def seq_reconstraction(query,pos_list):
            target_seq = "" 
            pre_pos    = min([pair[0] if pair[0] != None else 1000 for pair in pos_list])
            cigalet    = ""
            for i, pair in enumerate(pos_list):
                if pair[0] == None:
                    target_seq += "-"
                elif pair[0] - pre_pos > 1:
                    cigalet    += str(i-1) + query[pre_pos+1:pair[0]] 
                    target_seq += query[pair[0]] 
                    pre_pos    =  pair[0]
                else:
                    target_seq += query[pair[0]] 
                    pre_pos    = pair[0]
            return target_seq, cigalet 
                
        if method == "Needle":
            identity_list = []
            fastq_dict = read_fastq("{}.fastq".format(target))
            with open("{}_{}_glaln.csv".format(target,fname),"w") as o:
                count    = 0
                pair_seq = [] 
                pair_id  = [] 
                fasta_dict = collections.OrderedDict()
                for line in open("{}.needle".format(target)):
                    if line[0] == ">":
                        if count > 0:
                            pair_id.append(seqid[1:])
                            pair_seq.append(fasta_dict[seqid]) 

                        if count % 2 == 0 and count > 0:
                            #print(pair_seq[0])
                            #print(pair_seq[1])
                            #print("") 
                            identity = (len(pair_seq[0]) - Levenshtein.hamming(pair_seq[0],pair_seq[1])) * 1.0 / len(pair_seq[0].replace("-",""))
                            identity_list.append(identity * 100)
                            #print pair_seq[0]
                            #print pair_seq[1] 
                            #print (len(pair_seq[0]) - Levenshtein.hamming(pair_seq[0],pair_seq[1])) * 1.0 / len(pair_seq[0].replace("-",""))
                            s = 0
                            reference = pair_seq[0] 
                            if "-" not in reference and s > -1:
                                region = pair_seq[1]
                                o.write(",".join([seqid[1:],"",region]) + "\n")
                            else:
                                s = pair_seq[0].replace("-","").find(self.regions["sequences"][0]) 
                                reference = pair_seq[0].replace("-","") 
                                if s > -1:
                                    for n in range(len(pair_seq[0])):
                                        s = pair_seq[0][n:].replace("-","").find(reference) 
                                        if s == 0:
                                            break        
                                    s = n

                                    for n in range(s,len(pair_seq[0])):
                                        m = pair_seq[0][:n].replace("-","").find(reference) 
                                        if m > -1:
                                            break
                                    e = n
                                    _ins = ""
                                    char_count = 0 
                                    query   = pair_seq[1][s:e]
                                    asciies = fastq_dict[pair_id[1]][3][s:e]
                                    #print(len(query),len(asciies),asciies) 
                                    scores = [] 
                                    for asc in asciies:
                                        score = ord(asc)-33
                                        scores.append(score) 
                                    #print(scores) 
                                    region_scores = []
                                    region = ""
                                    flag_r   = 0 
                                    flag_q   = 0 
                                    #num += 1 
                                    pos = 0
                                    #print(query) 
                                    #print(pair_seq[0][s:e]) 
                                    #print(len(scores))
                                    for q, r in zip(query,pair_seq[0][s:e]):
                                        if r != "-":
                                            flag_r = 1
                                        
                                        if r == "-" and flag_r == 1:
                                            _ins += (str(char_count) + q)
                                            region_scores.append(scores[pos]) 

                                        elif r != "-":
                                            region += q
                                            char_count += 1
                                            #Eprint(len(scores))
                                            if q != "-":
                                                region_scores.append(scores[pos]) 
                                        
                                        if q == "-":
                                            pass  
                                        else:
                                            pos += 1
                                    #print(region_scores) 
                                    o.write(",".join([seqid[1:],_ins,region,str(identity),str(min(region_scores)),str(sum(region_scores)*1.0/len(region_scores))]) + "\n")
                                        
                                else:
                                    pass 
                            pair_seq = [] 
                            pair_id  = []                       
                        count += 1
                        seqid = line.rstrip().replace(":","_").split(" ")[0] 
                        fasta_dict[seqid] = ""

                    elif line[0] == "#" or line[0] == "\n":
                        pass 
                    
                    else:
                        fasta_dict[seqid] += line.rstrip() 
           
            fig = plt.figure(figsize=(4,3))
            ax = fig.add_axes([0.1,0.1,0.8,0.8])
            a, b = np.histogram(identity_list,bins=100,range=(0,100),density=True)
            ax.hist(identity_list,bins=100,range=(0,100),color="#303030")
            ax.tick_params(labelsize=12,pad=4)
            #ax.set_ylim(0,1)
            ax.set_xlim(0,100)
            ax.set_xlabel("Identity(%)",fontsize=12)
            ax.set_ylabel("#Read",fontsize=12)
            ax2 = ax.twinx()
            cumsum = [0] + list(np.cumsum(a))
            ax2.tick_params(labelsize=12,pad=4)
            ax2.plot(cumsum,color="#AAAAAA")      
            ax2.set_xlim(0,100)
            ax2.set_ylim(0,1) 
            ax2.set_ylabel("Cumulative ratio",fontsize=12)
            fig.savefig("Read_info.pdf",bbox_inches="tight") 
        else:
            pass 

    def barcode_count(self, target, visualization=True, lim_quality=0):
        min_match = self.min_match
        barcode_dict  = collections.defaultdict(dict) 
        reference_seq = self.references[0]
        for i, barcode_info in enumerate(self.targets):
            barcode_info = tuple(list(map(lambda x : tuple(x) if type(x) == list else x ,barcode_info)))
            with open("{}_target_region_glaln.csv".format(target)) as f: 
                barcode_ref = reference_seq[barcode_info[2][0]:barcode_info[2][1]]
                barcode_dict[barcode_info] = collections.defaultdict(int)
                with open("barcode{}_list.csv".format(i),"w") as o:
                    for line in f:
                        _id, ins, query, identity, min_quality, avg_quality = line.rstrip().split(",")
                        identity = float(identity) 
                        min_quality = int(min_quality)
                        if identity > min_match and min_quality > lim_quality:
                            barcode_query = query[barcode_info[2][0]:barcode_info[2][1]]
                            if len(ins) > 0:
                                pattern1  = r"([0-9]*)"
                                pattern2  = r"([^0-9]*)"
                                positions = re.findall(pattern1,ins)[:-2][0::2]  
                                ins_pos   = list(map(int,positions))
                                ins_char  = list(re.findall(pattern2,ins)[:-1][2::3])        
                                ins_num   = 0
                                for pos, char in zip(ins_pos,ins_char):
                                    pos = pos - barcode_info[2][0]
                                    if 0 < pos < (barcode_info[2][1] - barcode_info[2][0]):
                                        barcode_query = barcode_query[0:pos+ins_num] + char +  barcode_query[pos+ins_num:]
                                        ins_num += 1
                                    else:
                                        pass 
                            else:
                                pass 
                            o.write("{},{}\n".format(_id,barcode_query)) 
                            barcode_dict[barcode_info][barcode_query] += 1
            
        for i, barode_info in enumerate(barcode_dict):
            with open("barcode{}_count.csv".format(i),"w") as o:
                for barcode_seq in barcode_dict[barcode_info]:
                    o.write("{},{}\n".format(barcode_seq,barcode_dict[barcode_info][barcode_seq]))

        return barcode_dict
    
    def calc_editing_spectrum(self, target, visualization=True, be=False, lim_quality=0):
        min_match = self.min_match
        def colorbar(ax,color_dict,ref_seq,char=False):
            bars = ax.bar(list(range(len(ref_seq))), [1] * (len(ref_seq)), width=1.0, edgecolor="#BBBBBB", linewidth=0.2, align="edge")
            ax.set_xlim(0,len(ref_seq))
            ax.set_ylim(0,1)
            p = 0
            for bar, c in zip(bars,ref_seq):
                color = color_dict[c]
                bar.set_facecolor(color)
                if char == True:
                    ax.text(p+0.5,0.4,c,va="center",ha="center",fontsize=9,zorder=100)     
                p += 1 
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.patch.set_alpha(0.0)
        
        def plot_editing_ratio(ax,ratio,ymax):
            ax.bar(np.array(list(range(len(ratio)))) + 0.5, ratio, lw=0.0, width=0.9, facecolor="#303030", edgecolor="#999999", align="center") 
            ax.tick_params(labelsize=14, pad=2)
            ax.set_xlim(0,len(ratio))
            ax.set_ylim(0,ymax)
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            #ax.set_title(sample_dict[key1]["base_editor"], fontsize=16)
        
        def plot_targets(positions, reference_seq, ranges):
            poscolors_FW = ["#FFF2F2"] * len(reference_seq) 
            poscolors_RV = ["#F2F2FF"] * len(reference_seq)
            for target_info in targets:
                target_positions_list_FW = []
                target_positions_list_RV = []
                s = target_info[2][0] - ranges[0][0]
                e = target_info[2][1] - ranges[0][0]
                if target_info[1] == 1:
                    target_positions_list_FW.append(list(range(s,e))) 
                else:
                    target_positions_list_RV.append(list(range(s,e))) 

                for p in range(len(reference_seq)):
                    if sum([1 if p in pos else 0 for pos in target_positions_list_FW]) > 0:
                        poscolors_FW[p] = target_info[3]
                    else:
                        pass
                 
                for p in range(len(reference_seq)):
                    if sum([1 if p in pos else 0 for pos in target_positions_list_RV]) > 0:
                        poscolors_RV[p] = target_info[3]
                    else:
                        pass

            ax  = fig.add_axes(positions) 
            ax.spines["right"].set_visible(False)
            #ax.spines["bottom"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)
                        
            ax.tick_params(pad=2, labelsize=14)
            bars = ax.bar(list(range(len(reference_seq))), [1] * (len(reference_seq)), bottom=1,  width=1.0, linewidth=0.0, align="edge")
            for bar, color in zip(bars,poscolors_FW):
                bar.set_facecolor(color)
                bar.set_edgecolor(color)
           
            bars = ax.bar(list(range(len(reference_seq))), [1] * (len(reference_seq)), bottom=0,  width=1.0, linewidth=0.0, align="edge")
            for bar, color in zip(bars,poscolors_RV):
                bar.set_facecolor(color)
                bar.set_edgecolor(color)
            
            ax.set_xlim(0,len(poscolors_FW))
            if self.position_type == "relative":
                pass
            else:
                locs = ax.xaxis.get_ticklocs()   
                ax.set_xticklabels(list(map(lambda x:str(int(x+self.regions["ranges"][0][0] - int(self.position_type))) ,locs)))

            ax.set_xlabel("Position (bp)",fontsize=14)
            ax.set_ylim(0,2)
            ax.set_yticks([])
            return poscolors_FW, poscolors_RV

        def plot_pattern(ax1, ax2, ax3, pattern_dict, targets, ranges, control, view_num=15, count=0, length=None):
            pos_sets = []
            for target_info in targets:
                pos_info = target_info[2] 
                pos_sets.append(pos_info)  
            pos_sets = list(zip(*pos_sets)) 
            s,e = min(pos_sets[0])-ranges[0][0], max(pos_sets[1])-ranges[0][0]
            
            pattern_combi = list(pattern_dict.items())
            pattern_combi.sort(key=lambda x:x[1])
            pattern_combi.reverse()
            m = view_num
            new_pattern_combi = [] 
            
            for key, value in pattern_combi:
                #print key
                if 1 in key or 2 in key:
                    new_pattern_combi.append([key,value]) 
                else:
                    pass 
            
            ins_count    = 0
            del_count   = 0
            indel_count = 0 
            for key, value in new_pattern_combi:
                if 1 in key[s:e]:
                    del_count += value
                if 2 in key[s:e]:
                    ins_count += value
                if 1 in key[s:e] or 2 in key[s:e]: 
                    indel_count += value
            
            ins_count   = ins_count * 1.0 / count
            del_count   = del_count * 1.0 / count
            indel_count = indel_count* 1.0 / count 
            
            with open("indel_summary.csv","w") as o:
                o.write("indel_ratio,del_ratio,ins_ratio\n") 
                o.write(",".join(list(map(str,[indel_count,del_count,ins_count]))) + "\n")
            fig_summary = plt.figure(figsize=(3,3))
            ax_summary  = fig_summary.add_axes([0.1,0.1,0.8,0.8])
            ax_summary.tick_params(labelsize=16)
            ax_summary.bar([0,1,2],[indel_count,del_count,ins_count],width=0.8) 
            ax_summary.set_xlim(-0.5,2.5) 
            ax_summary.set_ylim(0,math.ceil(max([indel_count,del_count,ins_count]) * 10.0) * 1.0 / 10) 
            ax_summary.set_xticks([0,1,2]) 
            ax_summary.set_xticklabels(["indel ratio","del ratio","ins ratio"],rotation=45,ha="right")
            ax_summary.spines["top"].set_visible(False) 
            ax_summary.spines["right"].set_visible(False) 
            fig_summary.savefig("indel_summary.pdf",bbox_inches="tight")
            
            if control[1] != "control" and control[1] != "N.A.":
                with open("../../{}/{}/indel_summary.csv".format(control[0],control[1])) as f:
                    f.readline() 
                    con_ratios = list(map(float,f.readline().split(",")))
                fig_summary = plt.figure(figsize=(3,3))
                ax_summary  = fig_summary.add_axes([0.1,0.1,0.8,0.8])
                ax_summary.tick_params(labelsize=16)
                ax_summary.bar([0,1,2],[indel_count-con_ratios[0],del_count-con_ratios[1],ins_count-con_ratios[2]],width=0.8) 
                with open("indel_summary_normalized.csv","w") as o:
                    o.write("indel_ratio,del_ratio,ins_ratio\n") 
                    o.write(",".join(list(map(str,[indel_count-con_ratios[0],del_count-con_ratios[1],ins_count-con_ratios[2]]))) + "\n")
                    
                ax_summary.set_xlim(-0.5,2.5) 
                ax_summary.set_ylim(0,math.ceil(max([indel_count-con_ratios[0],del_count--con_ratios[1],ins_count-con_ratios[2]]) * 10.0) * 1.0 / 10) 
                ax_summary.set_xticks([0,1,2]) 
                ax_summary.set_xticklabels(["indel ratio","del ratio","ins ratio"],rotation=45,ha="right")
                ax_summary.spines["top"].set_visible(False) 
                ax_summary.spines["right"].set_visible(False) 
                fig_summary.savefig("indel_summary_normalized.pdf",bbox_inches="tight")


            for key, value in new_pattern_combi[0:view_num]:
                n = 0 
                for c in key:
                    if c == 1:
                        facecolor = "#FF0000"
                    elif c == 2:
                        facecolor = "#0000FF"
                    else:
                        facecolor = "#EEEEEE"
                    ax1.barh([m],[1.0], left=n, lw=0.0, height=0.9, align="center", facecolor=facecolor) 
                    n += 1
                ax2.barh([m], [value], lw=0.0, height=0.9, align="center",  facecolor="#303030") 
                m -= 1
            
            ax1.spines["right"].set_visible(False)
            ax1.spines["bottom"].set_visible(False)
            ax1.spines["left"].set_visible(False)
            ax1.spines["top"].set_visible(False)
            ax1.set_ylim(0.5, view_num+0.5)
            ax1.set_xlim(0.0, length) 
            ax1.set_xticks([]) 
            ax1.set_yticks([]) 
            ax1.text(-0.02,0.5,"Indel pattern",fontsize=16,rotation=90,transform=ax1.transAxes,ha="center",va="center")
            
            ax2.tick_params(labelsize=16, pad=-18, direction="in") 
            ax2.spines["right"].set_visible(False)
            ax2.spines["top"].set_visible(False)
            ax2.spines["left"].set_visible(False)
            ax2.set_ylim(0.5, view_num+0.5)
            ax2.set_yticks([]) 
            ax2.spines["bottom"].set_position(("axes",1.02))
            ax2.text(0.5,1.25,"#Reads",fontsize=16,va="center",ha="center",transform = ax2.transAxes)

            m = view_num
            for key, value in new_pattern_combi[0:view_num]:
                if count == 0:
                    pass 
                else:
                    value = value * 1.0 / count
                ax3.barh([m], [value], lw=0.0, height=0.9,  align="center", facecolor="#303030", alpha=0.0)
            
            ax3.tick_params(labelsize=16, pad=3)  
            ax3.spines["right"].set_visible(False)
            ax3.spines["top"].set_visible(False)
            ax3.spines["left"].set_visible(False)
            ax3.set_ylim(0.5, view_num+0.5)
            ax3.set_yticks([]) 
            ax3.patch.set_alpha(0.0) 
            ax3.spines["bottom"].set_position(("axes",-0.02))
            ax3.set_xlabel("Frequency",fontsize=16)
            labels = list(map(str,ax3.get_xticks()))
            labels[0] = "0"
            ax3.xaxis.set_ticklabels(labels)
    
        def plot_substitution_pattern(ax1, ax2, ax3, ref_seq, pattern_dict, color_dict, control, view_num=8, count=0, deletion=True, length=None):
            pattern_combi = list(pattern_dict.items())
            pattern_combi.sort(key=lambda x:x[1])
            pattern_combi.reverse()
            m = view_num
                       
            if deletion == False:
                new_pattern_combi = [] 
                for key, value in pattern_combi:
                    if "-" not in key:
                        new_pattern_combi.append([key,value]) 
                    else:
                        pass  
                pattern_combi = new_pattern_combi
                #print(pattern_combi)
            
            new_pattern_combi = [] 
            for key, value in pattern_combi:
                if key != ref_seq:
                    new_pattern_combi.append([key,value]) 
                else:
                    pass
            
            pattern_combi = new_pattern_combi
            subpat_dict = collections.OrderedDict() 
            for c1 in  ["A","T","G","C"]:
                for c2 in ["A","T","G","C"]:
                    if c1 == c2:
                        pass
                    else:
                        subpat_dict[(c1,c2)] = 0 
            subpat_dict[("all",)] = 0  
            for key, value in pattern_combi:
                mut_set = []
                flag = 0
                for q, r in list(set(list(zip(key,ref_seq)))):
                    mut_set.append((q,r)) 
                    if q != r and q != "-" and q != "N":
                        subpat_dict[(r,q)] += value
                        flag = 1
                if flag == 1:
                    subpat_dict[("all",)] += value
            
            for key in subpat_dict:
                subpat_dict[key] = subpat_dict[key] * 1.0 / count                    
            
            with open("substitution_summary.csv","w") as o:
                o.write(",".join(list(map(lambda key:" to ".join(key), subpat_dict.keys()))) + "\n") 
                o.write(",".join(list(map(str,subpat_dict.values()))) + "\n")
            
            fig_summary = plt.figure(figsize=(8,3))
            ax_summary  = fig_summary.add_axes([0.1,0.1,0.8,0.8])
            ax_summary.tick_params(labelsize=16)
            ax_summary.bar(range(13),subpat_dict.values(),width=0.7) 
            ax_summary.set_xlim(-0.5,12.5) 
            ax_summary.set_ylim(0,math.ceil(max(subpat_dict.values()) * 10.0) * 1.0 / 10) 
            ax_summary.set_xticks(range(13)) 
            ax_summary.set_xticklabels(list(map(lambda key:" to ".join(key), subpat_dict.keys())),rotation=45,ha="right")
            ax_summary.spines["top"].set_visible(False) 
            ax_summary.spines["right"].set_visible(False) 
            fig_summary.savefig("substitution_summary.pdf",bbox_inches="tight")
            
            if control[1] != "control" and control[1] != "N.A.":
                with open("../../{}/{}/substitution_summary.csv".format(control[0],control[1])) as f:
                    f.readline() 
                    con_ratios = list(map(float,f.readline().split(",")))
                nor_ratios = list(map(lambda x,y:x-y, subpat_dict.values(), con_ratios))
                with open("substitution_summary_normalized.csv","w") as o:
                    o.write(",".join(list(map(lambda key:" to ".join(key), subpat_dict.keys()))) + "\n") 
                    o.write(",".join(list(map(str,nor_ratios))) + "\n")
                fig_summary = plt.figure(figsize=(8,3))
                ax_summary  = fig_summary.add_axes([0.1,0.1,0.8,0.8])
                ax_summary.tick_params(labelsize=16)
                ax_summary.bar(range(13),nor_ratios,width=0.7) 
                ax_summary.set_xlim(-0.5,12.5) 
                ax_summary.set_ylim(0,max([0.1,math.ceil(max(nor_ratios) * 10.0) * 1.0 / 10])) 
                ax_summary.set_xticks(range(13)) 
                ax_summary.set_xticklabels(list(map(lambda key:" to ".join(key), subpat_dict.keys())),rotation=45,ha="right")
                ax_summary.spines["top"].set_visible(False) 
                ax_summary.spines["right"].set_visible(False) 
                fig_summary.savefig("substitution_summary_normalized.pdf",bbox_inches="tight")
                 
            for key, value in pattern_combi[0:view_num]:
                n = 0 
                for q, r in zip(key,ref_seq):
                    if q == r:
                        facecolor = color_dict[r]
                        ax1.barh([m],[1.0], left=n, lw=0.5, height=0.9, align="center", edgecolor="#BBBBBB",facecolor=facecolor,zorder=10) 
                        ax1.text(n+0.5,m-0.1,r,va="center",ha="center",fontsize=9,zorder=10)     
                    else:
                        facecolor = color_dict[q]
                        ax1.barh([m],[1.0], left=n, lw=0.5, height=0.9, align="center", edgecolor="k", facecolor=facecolor,zorder=100) 
                        ax1.text(n+0.5,m-0.1,q,fontweight="bold",va="center",ha="center",fontsize=9,zorder=100) 
                    n += 1
                ax2.barh([m], [value], lw=0.0, height=0.9, align="center",  facecolor="#303030") 
                m -= 1

            ax1.spines["right"].set_visible(False)
            ax1.spines["bottom"].set_visible(False)
            ax1.spines["left"].set_visible(False)
            ax1.spines["top"].set_visible(False)
            ax1.set_ylim(0.5, view_num+0.5)
            ax1.set_xlim(0.0, length) 
            ax1.set_xticks([]) 
            ax1.set_yticks([]) 
            ax1.text(-0.04*23.0/len(ref_seq),0.55,"Read pattern",fontsize=14,rotation=90,transform=ax1.transAxes,ha="center",va="center")
            
            ax2.tick_params(labelsize=14, pad=-14, direction="in") 
            ax2.spines["right"].set_visible(False)
            ax2.spines["top"].set_visible(False)
            ax2.spines["left"].set_visible(False)
            ax2.set_ylim(0.5, view_num+0.5)
            ax2.set_yticks([]) 
            ax2.spines["bottom"].set_position(("axes",1.02))
            ax2.text(0.5,1.27,"#Reads",fontsize=14,va="center",ha="center",transform = ax2.transAxes)
            
            m = view_num
            for key, value in pattern_combi[0:view_num]:
                if count == 0:
                    pass 
                else:
                    value = value * 1.0 / count
                ax3.barh([m], [value], lw=0.0, height=0.9, align="center",  facecolor="#303030") 
                m -= 1
            
            ax3.tick_params(labelsize=14, pad=0)  
            ax3.spines["right"].set_visible(False)
            ax3.spines["top"].set_visible(False)
            ax3.spines["left"].set_visible(False)
            ax3.set_ylim(0.5, view_num+0.5)
            ax3.set_yticks([]) 
            ax3.patch.set_alpha(0.0) 
            ax3.spines["bottom"].set_position(("axes",-0.02))
            labels = list(map(str,ax3.get_xticks()))
            labels[0] = "0.0"
            ax3.xaxis.set_ticklabels(labels)
            ax3.set_xlabel("Frequency",fontsize=14)

        s = self.regions["ranges"][0][0]
        left_margin   = self.regions["ranges"][0][0]
        reference_seq = self.references[0][left_margin:self.regions["ranges"][0][1]]  
        targets       = self.targets
        ins_rate = [0] * len(reference_seq)
        del_rate = [0] * len(reference_seq)
        mut_rate = [0] * len(reference_seq) 
        del_pattern_dict = collections.defaultdict(int)
        ins_pattern_dict = collections.defaultdict(int)
        del_ins_pattern_dict = collections.defaultdict(int)
        nuc_occupancy = [{"A":0,"T":0,"G":0,"C":0,"N":0,"-":0} for i in range(len(reference_seq))]  
         
        if be == True:
            mut_pattern_dict = collections.OrderedDict()
            for target_info in targets:
                if type(target_info[3]) == list:
                    target_info[3] = tuple(target_info[3])
                target_info[2] = tuple(target_info[2])
                mut_pattern_dict[tuple(target_info)] = collections.defaultdict(int)
        count = 1
        for line in open("{}_target_region_glaln.csv".format(target)):
            _id, ins, query, identity, min_quality, avg_quality = line.rstrip().split(",")
            identity = float(identity) 
            min_quality = int(min_quality) 
            if identity > min_match and min_quality > lim_quality:
                query = query[self.regions["ranges"][0][0]:self.regions["ranges"][0][1]]
                del_pattern = query.translate(str.maketrans("ATGCN-","000001")) 
                del_pattern = list(map(int,del_pattern))
                del_pattern_dict[tuple(del_pattern)] += 1     
                ins_pattern = [0] * len(reference_seq)
                #del_ins_pattern_dict[tuple(del_pattern)] += 1
                
                if len(ins) > 0:
                    pattern1  = r"([0-9]*)"
                    pattern2  = r"([^0-9]*)"
                    positions = re.findall(pattern1,ins)[:-2][0::2]  
                    ins_pos = set(list(map(int,positions)))
                    for pos in ins_pos:
                        if 0 < pos - self.regions["ranges"][0][0] < len(reference_seq):
                            ins_pattern[pos - self.regions["ranges"][0][0]] = 1 
                            ins_rate[pos - self.regions["ranges"][0][0]] += 1
                    ins_pattern_dict[tuple(ins_pattern)] += 1
                
                if be == True:
                    for target_info in mut_pattern_dict:
                        s = target_info[2][0] - self.regions["ranges"][0][0]
                        e = target_info[2][1] - self.regions["ranges"][0][0]
                        mut_pattern_dict[tuple(target_info)][query[s:e]] += 1 
                
                del_ins_pattern = [1 if de == 1 else 2 if ins == 1 else 0 for de, ins in zip(del_pattern, ins_pattern)]
                del_ins_pattern_dict[tuple(del_ins_pattern)] += 1
                #for pos, nuc in zip(re.findall(pattern1,ins),re.findall(pattern2,ins)):
                pos = 0
                for q, r in zip(query, reference_seq):
                    if q == r:
                        pass 
                    else:
                        if q == "-":
                            del_rate[pos] += 1
                        else:
                            mut_rate[pos] += 1
                            nuc_occupancy[pos][q] += 1
                    pos += 1
                count += 1
        print("Read count:{}".format(count))
        ins_rate = list(map(lambda x: x*1.0/count, ins_rate))
        del_rate = list(map(lambda x: x*1.0/count, del_rate))
        mut_rate = list(map(lambda x: x*1.0/count, mut_rate))
        self[target]["count"] = count 
        #self[target]["ins_rate"] = ins_rate
        #self[target]["del_rate"] = del_rate
        #self[target]["mut_rate"] = mut_rate
        
        #print(self.regions["ranges"][0])  
        if self.control[1] != "control" and self.control[1] != "N.A.":
            with open("../../{}/{}/editing_spectrums.json".format(self.control[0],self.control[1])) as f:
                control_spectrums = json.load(f) 
                self.conobj[target]["count"] = control_spectrums["sample"]["count"] 
                self.conobj[target]["mut_rate"] = control_spectrums["sample"]["mut_rate"] 
                self.conobj[target]["ins_rate"] = control_spectrums["sample"]["ins_rate"]
                self.conobj[target]["del_rate"] = control_spectrums["sample"]["del_rate"]

            with open("editing_spectrums.json","w") as o:
                editing_spectrums = {"sample":{},"control":{},"normalized":{}}
                editing_spectrums["sample"]["count"] = count
                editing_spectrums["sample"]["mut_rate"] = mut_rate
                editing_spectrums["sample"]["ins_rate"] = ins_rate
                editing_spectrums["sample"]["del_rate"] = del_rate
                editing_spectrums["sample"]["nuc_occupancy"] = nuc_occupancy         
                editing_spectrums["control"]["count"] = self.conobj[target]["count"]
                editing_spectrums["control"]["mut_rate"] = self.conobj[target]["mut_rate"]
                editing_spectrums["control"]["ins_rate"] = self.conobj[target]["ins_rate"]
                editing_spectrums["control"]["del_rate"] = self.conobj[target]["del_rate"]
                #print(self.conobj.__dict__) 
                #ins_rate = list(map(lambda x,y: x-y if x-y > 0 else 0.0, ins_rate, self.conobj[target]["ins_rate"][self.regions["ranges"][0][0]:self.regions["ranges"][0][1]]))
                #del_rate = list(map(lambda x,y: x-y if x-y > 0 else 0.0, del_rate, self.conobj[target]["del_rate"][self.regions["ranges"][0][0]:self.regions["ranges"][0][1]])) 
                #mut_rate = list(map(lambda x,y: x-y if x-y > 0 else 0.0, mut_rate, self.conobj[target]["mut_rate"][self.regions["ranges"][0][0]:self.regions["ranges"][0][1]]))  
                print(list(zip(mut_rate,self.conobj[target]["mut_rate"]))) 
                ins_rate = list(map(lambda x,y: x-y if x-y > 0 else 0.0, ins_rate, self.conobj[target]["ins_rate"]))
                del_rate = list(map(lambda x,y: x-y if x-y > 0 else 0.0, del_rate, self.conobj[target]["del_rate"])) 
                mut_rate = list(map(lambda x,y: x-y if x-y > 0 else 0.0, mut_rate, self.conobj[target]["mut_rate"])) 
                print(mut_rate)
                editing_spectrums["normalized"]["mut_rate"] = mut_rate
                editing_spectrums["normalized"]["ins_rate"] = ins_rate
                editing_spectrums["normalized"]["del_rate"] = del_rate
                json.dump(editing_spectrums,o) 

        else:
            with open("editing_spectrums.json","w") as o:
                editing_spectrums = {"sample":{},"control":{},"normalized":{}}
                editing_spectrums["sample"]["count"] = count
                editing_spectrums["sample"]["mut_rate"] = mut_rate
                editing_spectrums["sample"]["ins_rate"] = ins_rate
                editing_spectrums["sample"]["del_rate"] = del_rate
                editing_spectrums["sample"]["nuc_occupancy"] = nuc_occupancy         
                json.dump(editing_spectrums,o)

        #print(ins_rate,del_rate,mut_rate)
        max_rate = max(list(map(max,[ins_rate,del_rate,mut_rate])))
        max_rate = max_rate if max_rate > 0 else 0.01
        max_rate = round(max_rate + 0.5*(pow(10,(int(math.log10(max_rate)) - 1))),abs(int(math.log10(max_rate)) - 1))
        max_rate = 0.5
        oc_dict = {"A":[],"T":[],"G":[],"C":[],"N":[],"-":[]}
        for j in range(len(nuc_occupancy)):
            for key in sorted(list(nuc_occupancy[j].keys())):
                oc_dict[key].append(nuc_occupancy[j][key]) 
        
        with open("indel_pattern.txt","w") as o:
            o.write("\t".join(["indel_pattern(1=deletion,2=insertion)","Read count","Frequency"]) + "\n") 
            o.write("\t".join([reference_seq,"N.A.","N.A."]) + "\n") 
            for pattern in list(reversed(list(sorted(list(del_ins_pattern_dict.keys()),key=lambda x:del_ins_pattern_dict[x])))):
                o.write("\t".join(["".join(list(map(str,pattern))),str(del_ins_pattern_dict[pattern]),str(del_ins_pattern_dict[pattern]*1.0/count)]) + "\n") 
        
        num = 1
        for target_info in mut_pattern_dict: 
            with open("mutation_pattern_target{}.txt".format(num),"w") as o:
                for pattern in list(reversed(list(sorted(list(mut_pattern_dict[target_info].keys()),key=lambda x:mut_pattern_dict[target_info][x])))):
                    o.write("\t".join([pattern,str(mut_pattern_dict[target_info][pattern])]) + "\n") 
            num += 1 

        if visualization == True:
            color_dict  = {"G":"#f2f059", "C":"#74b2d7", "A":"#79E5B7", "T":"#ff776c", "N":"#FFFFFF", "-":"#FFFFFF"}
            #Each editing spectrum         
            fig = plt.figure(figsize=(4.5*(self.regions["ranges"][0][1]-self.regions["ranges"][0][0])/60.0,2.4)) 
            ax  = fig.add_axes([0.1,1.1,0.8,0.5]) 
            #max_rate = max(mut_rate)
            #max_rate = max_rate if max_rate > 0 else 0.01
            #max_rate = round(max_rate + 0.5*(pow(10,(int(math.log10(max_rate)) - 1))),abs(int(math.log10(max_rate)) - 1))
            plot_editing_ratio(ax,mut_rate,max_rate)
            ax.set_ylabel("Frequency",fontsize=14)
            ax.set_title("Substitution spectrum",fontsize=14)
            if self.position_type == "relative":
                pass 
            else:
                locs = ax.xaxis.get_ticklocs()   
                ax.set_xticklabels(list(map(lambda x:str(int(x+self.regions["ranges"][0][0] - int(self.position_type))), locs)))
            ax.set_xticks([])
            ax.set_xticklabels([])
            poscolors_FW, poscolors_RV = plot_targets([0.1,1.015,0.8,0.08], reference_seq, self.regions["ranges"])

            ax  = fig.add_axes([0.1,2.1,0.8,0.5]) 
            #max_rate = max(del_rate)
            #max_rate = max_rate if max_rate > 0 else 0.01
            #max_rate = round(max_rate + 0.5*(pow(10,(int(math.log10(max_rate)) - 1))),abs(int(math.log10(max_rate)) - 1))
            plot_editing_ratio(ax,del_rate,max_rate)
            ax.set_ylabel("Frequency",fontsize=14)
            ax.set_title("Deletion spectrum",fontsize=14)
            if self.position_type == "relative":
                pass 
            else:
                locs = ax.xaxis.get_ticklocs()   
                ax.set_xticklabels(list(map(lambda x:str(int(x+self.regions["ranges"][0][0] - int(self.position_type))), locs)))
            ax.set_xticks([])
            ax.set_xticklabels([])
            poscolors_FW, poscolors_RV = plot_targets([0.1,2.015,0.8,0.08], reference_seq, self.regions["ranges"])


            ax  = fig.add_axes([0.1,3.1,0.8,0.5]) 
            #max_rate = max(ins_rate)
            #max_rate = max_rate if max_rate > 0 else 0.01
            #max_rate = round(max_rate + 0.5*(pow(10,(int(math.log10(max_rate)) - 1))),abs(int(math.log10(max_rate)) - 1))
            plot_editing_ratio(ax,ins_rate,max_rate)
            ax.set_ylabel("Frequency",fontsize=14)
            ax.set_title("Insertion spectrum",fontsize=14)
            if self.position_type == "relative":
                pass 
            else:
                locs = ax.xaxis.get_ticklocs()   
                ax.set_xticklabels(list(map(lambda x:str(int(x+self.regions["ranges"][0][0] - int(self.position_type))) ,locs)))
            ax.set_xticks([])
            ax.set_xticklabels([])
            poscolors_FW, poscolors_RV = plot_targets([0.1,3.015,0.8,0.08], reference_seq, self.regions["ranges"])

            ax  = fig.add_axes([0.1,0.17,0.8,0.46])     
            for j in range(len(nuc_occupancy)):
                bottom = 0
                sum_value = sum(list(nuc_occupancy[j].values()))
                if sum_value == 0:
                    sum_value = 0.1
                for key in sorted(list(nuc_occupancy[j].keys())):
                    if mut_rate[j] >= 0.01:
                        ax.bar([j+0.5], [nuc_occupancy[j][key] * 1.0 / sum_value], lw=0.0, width=0.9, facecolor=color_dict[key], edgecolor="#999999", align="center", bottom=bottom)
                    bottom += nuc_occupancy[j][key] * 1.0 / sum_value

            ax.tick_params(labelsize=14, pad=2)
            ax.set_xticks([]) 
            ax.set_yticks([])
            ax.set_xlim(0,len(reference_seq))
            ax.set_ylim(0,1.0) 
            ax.set_yticks([0,0.5,1.0]) 
            ax.set_ylabel("Occupancy",fontsize=14)
            ax.set_title("Mutation composition",fontsize=14)
                
            ax  = fig.add_axes([0.1,0.10,0.8,0.06]) 
            ax.tick_params(labelsize=14, pad=2)
            colorbar(ax,color_dict, reference_seq)
            ax.tick_params(axis="x", length=0)
            poscolors_FW, poscolors_RV = plot_targets([0.1,0.00,0.8,0.1],reference_seq, self.regions["ranges"])
            fig.patch.set_alpha(0.0) 
            fig.savefig("editing_spectrums_{}.pdf".format(target), bbox_inches="tight")
            if be == True:
                pass 
            
            fig = plt.figure(figsize=(7,2))
            ax1 = fig.add_axes([0.10,0.1,0.8,0.8])   
            ax2 = fig.add_axes([0.92,0.1,0.24,0.8],label="a") 
            ax3 = fig.add_axes([0.92,0.1,0.24,0.8],label="b") 

            plot_pattern(ax1, ax2, ax3, del_ins_pattern_dict, self.targets, self.regions["ranges"], self.control, count=count, length=len(reference_seq),) 
            ax  = fig.add_axes([0.1,0.02,0.8,0.07]) 
            ax.tick_params(labelsize=16, pad=2)
            colorbar(ax, color_dict, reference_seq)
            ax.tick_params(axis="x", length=0) 
            
            ax  = fig.add_axes([0.1,-0.08,0.8,0.1]) 
            ax.spines["right"].set_visible(False)
            #ax.spines["bottom"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["top"].set_visible(False)            
            ax.tick_params(pad=2, labelsize=16)
            bars = ax.bar(list(range(len(poscolors_FW))), [1] * (len(poscolors_FW)), bottom=1,  width=1.0, linewidth=0.0, align="edge")
            for bar, color in zip(bars,poscolors_FW):
                bar.set_facecolor(color)
                bar.set_edgecolor(color)
           
            bars = ax.bar(list(range(len(poscolors_RV))), [1] * (len(poscolors_RV)), bottom=0,  width=1.0, linewidth=0.0, align="edge")
            for bar, color in zip(bars,poscolors_RV):
                bar.set_facecolor(color)
                bar.set_edgecolor(color)
            ax.set_xlim(0,len(poscolors_FW))
            ax.set_ylim(0,2)
            ax.set_yticks([])
            if self.position_type == "relative":
                pass 
            else:
                locs = ax.xaxis.get_ticklocs()   
                ax.set_xticklabels(list(map(lambda x:str(int(x+self.regions["ranges"][0][0] - int(self.position_type))) ,locs)))
            ax.set_xlabel("Position (bp)",fontsize=16)
            fig.patch.set_alpha(0.0) 
            fig.savefig("deletion_pattern_{}.pdf".format(target), bbox_inches="tight")
            plt.close() 

            if be == True:
                num = 1 
                for target_info in mut_pattern_dict:
                    if int(target_info[-1]) == 1:
                        fig = plt.figure(figsize=(3*len(target_info[0])/23.0,1.6))
                        ax1 = fig.add_axes([0.1,0.1,0.8,0.8])   
                        ax2 = fig.add_axes([0.90+0.02*23.0/len(target_info[0]),0.1,0.6*18.0/len(target_info[0]),0.8],label="a") 
                        ax3 = fig.add_axes([0.90+0.02*23.0/len(target_info[0]),0.1,0.6*18.0/len(target_info[0]),0.8],label="b") 
                        if target_info[1] == 1:
                            plot_substitution_pattern(ax1, ax2, ax3, target_info[0], mut_pattern_dict[target_info], color_dict, self.control, count=count, length=len(target_info[0])) 
                        else:
                            #print((ax1, ax2, ax3, target_info[0].translate(str.maketrans("ATGC","TACG"))[::-1], mut_pattern_dict[target_info], color_dict, count, len(target_info[0]))) 
                            plot_substitution_pattern(ax1, ax2, ax3, target_info[0].translate(str.maketrans("ATGC","TACG"))[::-1], mut_pattern_dict[target_info], color_dict, self.control, count=count, length=len(target_info[0])) 
                        ax  = fig.add_axes([0.1,0.92,0.8,0.1]) 
                        if target_info[1] == 1:
                            colorbar(ax, color_dict, target_info[0], char=True)
                        else:
                            colorbar(ax, color_dict, target_info[0].translate(str.maketrans("ATGC","TACG"))[::-1], char=True)
                        ax.tick_params(axis="x", length=0, labelsize=14) 
                        ax.text(0.5,1.8,"Reference sequence",fontsize=14,va="center",ha="center",transform=ax.transAxes) 
                        fig.patch.set_alpha(0.0) 
                        fig.savefig("substitution_pattern_{}_target{}.pdf".format(target,num), bbox_inches="tight")
                        plt.close()
                    num += 1

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-s","--sample",  type=str, default="None", help=" ")
    p.add_argument("-c","--control", type=str, default="None", help=" ")
    p.add_argument("-m","--min_identity", type=float, default=0.5, help=" ")
    p.add_argument("-r","--read_type", type=str, default="merged", help=" ")
    p.add_argument("-a","--analysis", type=str, default="base_editing", help=" ")
    p.add_argument("-q","--quality", type=int, default=0, help=" ")
    #p.add_argument("-q","--quality", type=int, default=0, help=" ")
    
    args   = p.parse_args()
    method = args.analysis
    sample = SAMPLE() 
    with open(args.sample + "/sample.json") as f:
        for key,value in json.load(f).items():
            sample[key] = value
    print("start")     
    if args.control != "None":
        pwd = os.getcwd()
        control = SAMPLE()     
        with open(args.control + "/sample.json") as f:
            for key,value in json.load(f).items():
                control[key] = value
        control.min_match = args.min_identity
        os.chdir(args.control) 
        if method == "base_editing" or "genome_editing":
            control.merge() 
            control.alignment(args.read_type) 
            control.extract_target_regions(args.read_type, method="Needle", fname="target_region", lim_quality=0)
            control.calc_editing_spectrum(args.read_type, visualization=True, be=True, lim_quality=args.quality)
        elif method == "barcode_count":
            control.barcode_count(args.read_type, visualization=True, be=True,lim_quality=args.quality)
        
        os.chdir(pwd) 
        sample.conobj = control
        os.chdir(args.sample) 
        sample.min_match = args.min_identity
        sample.control = [args.control.split("/")[0],args.control.split("/")[1]]
        if args.analysis == "base_editing" or "genome_editing":
            sample.extract_target_regions(args.read_type, method="Needle", fname="target_region", lim_quality=0)
            sample.calc_editing_spectrum(args.read_type, visualization=True, be=True, lim_quality=args.quality)
        elif args.analysis == "barcode_count":
            sample.barcode_count(args.read_type, visualization=True, lim_quality=args.quality)
    else:
        if (sample.control[1] == "N.A." or sample.control[1] == "control") == False:
            control = SAMPLE()    
            #print(sample.control)
            with open("./" + "/".join(sample.control) + "/sample.json") as f:
                for key,value in json.load(f).items():
                    control[key] = value
            control.min_match = args.min_identity
            sample.conobj = control
        os.chdir(args.sample) 
        sample.min_match = args.min_identity
        print("hoge") 
        if args.analysis == "base_editing" or args.analysis == "genome_editing":
            sample.merge()
            sample.alignment(args.read_type)
            sample.extract_target_regions(args.read_type, method="Needle", fname="target_region", lim_quality=-50)
            sample.calc_editing_spectrum(args.read_type, visualization=True, be=True,lim_quality=args.quality)
        elif args.analysis == "barcode_count":
            print("hgoe") 
            sample.barcode_count(args.read_type, visualization=True,lim_quality=args.quality)

