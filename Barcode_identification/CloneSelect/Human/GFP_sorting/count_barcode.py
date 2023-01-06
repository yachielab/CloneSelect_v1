import os 
import sys
import string
import itertools
import collections
from Bio import SeqIO 

before_bc = "CTACTAGAGGATCTATTTCCGGTGAATTCCCG"
targets   = ["112-V3-BC8","112-V3-BC9","112-V3-BC16","112-V3-BC19","112-V9-BC18","112-V9-BC19","112-V9-BC21","112-V10-BC14","112-V10-BC15","112-V10-BC17","112-V10-BC25","112-V10-BC29","112-V10-BC31","112-V10-BC33","112-V4-BC7","112-V4-BC8"] 

def hamming(seq1, seq2):
    return sum([ 1 for s1, s2 in zip(seq1, seq2) if s1 != s2])

def extract_barcode(file_name, before_bc, direction):
    for read in SeqIO.parse(file_name,"fastq"):
        if direction == "R1":
            seq = str(read.seq)
        else:
            seq = str(read.seq).translate(string.maketrans("ATGC","TACG"))[::-1]  
        
        pos = seq.find(before_bc)
        if pos > -1:
            bar_start = pos + len(before_bc) 
            bar_seq   = seq[bar_start:bar_start+25]
            bar_seqs_dict[read.id][direction] = (bar_seq,read.letter_annotations["phred_quality"][bar_start:bar_start+25]) 
        else:
            pass

def count_barcode(bar_seqs, bar_dict, target):
    for bar_seq in bar_seqs:
        distance_list = [] 
        for key in bar_dict.keys():
            dist = hamming(bar_seq[:14], bar_dict[key]["seq"][:14]) 
            distance_list.append([dist, key]) 
        distance_list.sort() 
        
        key     = distance_list[0][1]
        seq_len = len(bar_dict[key]["seq"])
        ew = ["GA" if nuc == "G" else nuc for nuc in bar_dict[key]["seq"][14:]] 
        combi = list(map(lambda x: "".join(x), list(itertools.product(*ew)))) 
        if distance_list[0][0] <=1 and bar_seq[14:seq_len] in combi:
            bar_dict[distance_list[0][1]]["count1"] += 1 
            if distance_list[0][0] <=1 and bar_seq[seq_len-3:seq_len] == "ATG":
                bar_dict[distance_list[0][1]]["count2"] += 1 
    #print bar_dict 
    return bar_dict 
        
if __name__ == "__main__":
    target   = sys.argv[1].split("/")[-1].split("_")[0] 
    #Extract barcode sequences from FASTQ input 
    bar_seqs_dict = collections.defaultdict(dict) 
    extract_barcode(sys.argv[1], before_bc, "R1")
    extract_barcode(sys.argv[2], before_bc, "R2")
    bar_seqs = [] 
    for key in bar_seqs_dict:
        if len(bar_seqs_dict[key]) == 2:
            if sum(bar_seqs_dict[key]["R1"][1]) >= sum(bar_seqs_dict[key]["R2"][1]):
                bar_seqs.append(bar_seqs_dict[key]["R1"][0]) 
            else:
                bar_seqs.append(bar_seqs_dict[key]["R2"][0]) 
    
    #print len(bar_seqs) 
    #Assign each barcode sequences to a specific barcode
    barcode_dict = {} 
    with open(sys.argv[3]) as f:
        for line in f:
            elements = line.rstrip().split(",")
            barcode_dict[elements[0]] = {} 
            barcode_dict[elements[0]]["seq"]    = elements[2]
            barcode_dict[elements[0]]["count1"] = 0
            barcode_dict[elements[0]]["count2"] = 0
            barcode_dict[elements[0]]["first_pop"] = float(int(elements[-1]))

    
    barcode_dict = count_barcode(bar_seqs, barcode_dict, target)
    others = set(list(barcode_dict.keys())) - set(targets)
    others = list(others)
    others.sort() 
    targets.sort()
    orders = targets + others 
   
    print(orders) 
    #print len(orders) 
    line1 = [target]
    line2 = [target]
    line3 = [target] 
    line4 = [target] 
    for key in orders:
        line1.append(str(barcode_dict[key]["count1"])) 
        line2.append(str(barcode_dict[key]["count2"])) 
        line3.append(str(barcode_dict[key]["count1"]/barcode_dict[key]["first_pop"])) 
        line4.append(str(barcode_dict[key]["count2"]/barcode_dict[key]["first_pop"])) 
    
    with open("barcode_counts_rep2.csv","aw") as o:
        o.write(",".join(line1) + "\n")
    
    with open("barcode_counts_rep2_only_atg.csv","aw") as o:
        o.write(",".join(line2) + "\n")

    with open("barcode_counts_rep2_normalized.csv","aw") as o:
        o.write(",".join(line3) + "\n")
    
    with open("barcode_counts_rep2_only_atg_normalized.csv","aw") as o:
        o.write(",".join(line4) + "\n")
