import os 
import re
import sys
import gzip
import time 
import json
import argparse
import Levenshtein
import collections
import subprocess 
import multiprocessing as mp 
import pickle
import string
import math
import numpy as np
from individual import *

class ANALYSIS(object):
    def __init__(self):
        self.fastq_list   = [] 
        self.fasta_list   = []
        self.index_dict   = {} 
        self.index_list   = []
        self.primers_dict = {} 
        self.index_primer_dict = collections.defaultdict(lambda:{"FW":[],"RV":[]})  
        self.experimental_dict = {}  
    
    def __getitem__(self, key):
        return self.__dict__[key]
    
    def __setitem__(self, key, value):
        self.__dict__[key] = value 

    def load_setting(self, samplefile, constfile):
        if "primers" not in os.listdir("."): 
            os.mkdir("./primers") 
        
        with open(samplefile) as _in:
            _in.readline()
            for line in _in:
                elements          = line.rstrip().split(",")
                project_id        = elements[0]
                sample_id         = elements[1]
                sample            = SAMPLE()
                sample.method     = elements[2]
                sample.id         = sample_id
                sample.FW_indices = elements[3].split(":")
                sample.RV_indices = elements[4].split(":")
                sample.FW_primer  = elements[5].upper() 
                sample.RV_primer  = elements[6].upper()
                sample.references = elements[7].upper().split(":")
                sample.targets    = elements[8].upper().split(":")
                if len(elements[9].split("|")) == 1:
                    sample.control    = (project_id,elements[9]) 
                else:
                    sample.control    = (elements[9].split("|")[0],elements[9].split("|")[1]) 
                sample.R1         = {} 
                sample.R2         = {} 
                sample.merged     = {} 
                whole             = sample.references[0]  
                for i, target in enumerate(sample.targets):
                    print(target) 
                    if set(target.split("|")[0]) <= set(["A","T","G","C","N","S","W","a","t","g","c","n","s","w"]):
                        target_seq = target.split("|")[0]
                        s = whole.find(target_seq)
                        e = s + len(target_seq)
                        direction = 1
                        if s == -1:
                            s = whole.find(target_seq.translate(str.maketrans("ATGC","TACG"))[::-1])
                            direction = -1
                            if s == -1:
                                print(elements)
                                raise Exception("Whole sequence of sample{} not contains target sequences".format(elements[0])) 
                        if len(target.split("|")) >= 3:
                            sample.targets[i] = [target_seq, direction, (s,s+len(target_seq)), target.split("|")[1], target.split("|")[2]]
                        elif len(target.split("|")) >= 2:
                            sample.targets[i] = [target_seq, direction, (s,s+len(target_seq)), target.split("|")[1], 1]
                        else:
                            sample.targets[i] = [target_seq, direction, (s,s+len(target_seq)), colors1[i], 1]

                    else:
                        s,e = target.split("|")[0:2]  
                        s = int(s) 
                        e = int(e) 
                        if s < e:
                            target_seq = whole[s:e] 
                            direction = 1
                        else:
                            target_seq = whole[e:s].translate(str.maketrans("ATGC","TACG"))[::-1]
                            direction = -1
                            s,e = e,s
                        
                        if len(target.split("|")) >= 3:
                            sample.targets[i] = [target_seq, direction, (s,e), target.split("|")[2], target.split("|")[3]]
                        elif len(target.split("|")) >= 2:
                            sample.targets[i] = [target_seq, direction, (s,e), target.split("|")[2], 1]
                        else:
                            sample.targets[i] = [target_seq, direction, (s,e), colors1[i], 1]


                if len(elements) == 10:
                    sample.regions  = {"ranges":[(0,len(whole))],"sequences":[whole[0:len(whole)]]}
                
                else:
                    sample.regions  = {"ranges":[],"sequences":[]}
                    for position_set in elements[10].split(":"):
                        s, e = list(map(int, position_set.split("|")))  
                        sample.regions["ranges"].append((s,e)) 
                        sample.regions["sequences"].append(whole[s:e])
                    sample.position_type = "relative"

                if len(elements) < 12:
                    sample.position_type = "relative"  
                else:
                    sample.position_type = elements[-1]
                self.index_list.append((int(sample.FW_indices[0]),int(sample.RV_indices[0])))
                self.index_primer_dict[(int(sample.FW_indices[0]),int(sample.RV_indices[0]))]["FW"].append((project_id,sample_id,sample.FW_primer)) 
                self.index_primer_dict[(int(sample.FW_indices[0]),int(sample.RV_indices[0]))]["RV"].append((project_id,sample_id,sample.RV_primer)) 
                self.experimental_dict[(project_id,sample_id)] = sample 
            #print(self.index_primer_dict) 
            #print(self.experimental_dict) 
        with open("./primers/primers.fasta","w") as o:
            const_dict = {}
            for line in open(constfile):
                if line[0] == ">":
                    key = line[1:].rstrip() 
                    const_dict[key] = ""
                else:
                    const_dict[key] += line.rstrip()
            
            for key in const_dict:
                if "Index" in key or "index" in key:
                    self.index_dict[key] = const_dict[key]
                
                if "primer" in key or "Primer" in key:
                    self.primers_dict[key] = const_dict[key]
                    o.write(">" + key.rstrip() + "\n")
                    o.write(const_dict[key].rstrip() + "\n") 
        os.chdir("./primers")   
        os.system("makeblastdb -in primers.fasta -out primers -dbtype nucl -parse_seqids") 
        os.chdir("../") 

    def load_data(self, miseq_data, split=500000):
        if "fragmented_fastq" not in os.listdir("."):
            os.mkdir("fragmented_fastq")
        
        os.chdir("fragmented_fastq")
        self.fastq_list = os.listdir("./")
        if len(self.fastq_list) <= 1:
            self.fastq_list = []
            for fastq in os.listdir("../" + miseq_data):
                i = 0
                n = 0
                flag = 0 
                o = open(fastq.replace(".fastq","_0.fastq"),"w") 
                self.fastq_list.append(fastq.replace(".fastq","_{}.fastq".format(0)))  
                for line in open("../" + miseq_data + "/" + fastq):
                    if flag == 1:
                        flag = 0 
                        o = open(fastq.replace(".fastq","_{}.fastq".format(int(n/split))),"w")
                        self.fastq_list.append(fastq.replace(".fastq","_{}.fastq".format(int(n/split))))  
                    i += 1
                    o.write(line)
                    if i % 4 == 0:
                        n += 1
                        if n % split == 0 and n > 0:
                            flag = 1 
                            o.close()
        os.chdir("../") 

    def fastq2fasta(self):
        if "fragmented_fasta" not in os.listdir("."):
            os.mkdir("fragmented_fasta")
        commands = [] 
        for fastq in self.fastq_list:
            commands.append("cat ./fragmented_fastq/{} | awk \'NR % 4 == 1 {{split($0, a, \" \"); gsub(\":\",\"_\",a[1]); print \">\" substr(a[1],2)}} NR % 4 == 2 {{print $0}}\' > ./fragmented_fasta/{}".format(fastq,fastq.replace("fastq","fasta")))
            self.fasta_list.append(fastq.replace("fastq","fasta"))
        
        for sub_commands in [commands[i:i+50] for i in range(0,len(commands),50)]:
            os.system("|".join(sub_commands)) 

    def search_primer(self, sge=False, num_threads=32):
        if num_threads == "MAX":
            num_threads = mp.cpu_count()
        header = "qseqid qseq sseqid sseq sstrand pident qlen length mismatch gapopen qstart qend sstart send gaps evalue bitscore btop" 
        header = header.split(" ")
        pwd = os.getcwd()
        def read_fasta(fasta_name):
            seq_dict = {} 
            with open(fasta_name) as f: 
                for line in f:
                    if line[0] == ">":
                        key = line[1:].rstrip() 
                        #key = "_".join(key.split("_")[:-1])
                        seq_dict[key] = "" 
                    else:
                        seq_dict[key] += line.rstrip() 
            return seq_dict 

        def read_fastq(fastq_name):
            seq_dict = {}  
            if fastq_name.split(".")[-1] == "gz":
                f = gzip.open(fastq_name.replace("'","").replace("\\",""), mode="rt", encoding='utf-8')
            else:
                f = open(fastq_name.replace("'","").replace("\\",""))  
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

        def read_blast(blastout, blast_dict, direction="R1"):
            with open(blastout) as f:
                for line in f:
                    line   = line.rstrip().split("\t")
                    qseqid = "_".join(line[0].split("_")) 
                    sseqid = line[2] 
                    blast_dict[qseqid][direction]["hit"][sseqid].append(dict(list(zip(header, line))))            
            return blast_dict 
        
        def search_index(index_dict,blast_dict):
            index_list = [index for index in index_dict.keys() if "Index_" in index] 
            index_blast_dict = {} 
            for target in blast_dict.keys():
                R1_hit_dict = blast_dict[target]["R1"]
                R2_hit_dict = blast_dict[target]["R2"]
                
                flag1  = 0
                index1 = "0"
                if "PS1.0-primer" in list(R1_hit_dict["hit"].keys()):
                    if R1_hit_dict["hit"]["PS1.0-primer"][0]["sstrand"]  == "plus":
                        seq         = R1_hit_dict["seq"] 
                        index_start = int(R1_hit_dict["hit"]["PS1.0-primer"][0]["qstart"]) - 10 - (int(R1_hit_dict["hit"]["PS1.0-primer"][0]["sstart"])-1) 
                        index_end   = int(R1_hit_dict["hit"]["PS1.0-primer"][0]["qstart"]) - 1 - (int(R1_hit_dict["hit"]["PS1.0-primer"][0]["sstart"])-1) 
                        index_seq   = seq[index_start:index_end] 
                        for index in index_list:
                            if Levenshtein.distance(index_dict[index],index_seq) <= 1:
                                index1    = index
                                flag1     = 1
                                break     
                
                if flag1 == 0 and "PS2.0-primer" in list(R1_hit_dict["hit"].keys()):
                    if R1_hit_dict["hit"]["PS2.0-primer"][0]["sstrand"]  == "plus":
                        seq         = R1_hit_dict["seq"] 
                        index_start = int(R1_hit_dict["hit"]["PS2.0-primer"][0]["qstart"]) - 10 - (int(R1_hit_dict["hit"]["PS2.0-primer"][0]["sstart"])-1) 
                        index_end   = int(R1_hit_dict["hit"]["PS2.0-primer"][0]["qstart"]) - 1 - (int(R1_hit_dict["hit"]["PS2.0-primer"][0]["sstart"])-1) 
                        index_seq   = seq[index_start:index_end] 
                        for index in index_list:
                            if Levenshtein.distance(index_dict[index],index_seq) <= 1:
                                index1    = index
                                flag1     = 1
                                break     

                flag2  = 0 
                index2 = "0"
                if "PS2.0-primer" in list(R2_hit_dict["hit"].keys()):
                    if R2_hit_dict["hit"]["PS2.0-primer"][0]["sstrand"]  == "plus":
                        seq         = R2_hit_dict["seq"] 
                        index_start = int(R2_hit_dict["hit"]["PS2.0-primer"][0]["qstart"]) - 10 - (int(R2_hit_dict["hit"]["PS2.0-primer"][0]["sstart"])-1)
                        index_end   = int(R2_hit_dict["hit"]["PS2.0-primer"][0]["qstart"]) - 1 - (int(R2_hit_dict["hit"]["PS2.0-primer"][0]["sstart"])-1)
                        index_seq   = seq[index_start:index_end] 

                        for index in index_list:
                            if Levenshtein.distance(index_dict[index],index_seq) <= 1:
                                index2    = index
                                flag2     = 1
                                break
                
                if flag2 == 0 and "PS1.0-primer" in list(R2_hit_dict["hit"].keys()):
                    if R2_hit_dict["hit"]["PS1.0-primer"][0]["sstrand"]  == "plus":
                        seq         = R2_hit_dict["seq"] 
                        index_start = int(R2_hit_dict["hit"]["PS1.0-primer"][0]["qstart"]) - 10 - (int(R2_hit_dict["hit"]["PS1.0-primer"][0]["sstart"])-1) 
                        index_end   = int(R2_hit_dict["hit"]["PS1.0-primer"][0]["qstart"]) - 1 - (int(R2_hit_dict["hit"]["PS1.0-primer"][0]["sstart"])-1) 
                        index_seq   = seq[index_start:index_end] 
                        for index in index_list:
                            if Levenshtein.distance(index_dict[index],index_seq) <= 1:
                                index2    = index
                                flag1     = 1
                                break     

                if (flag1 == 1 or flag2 == 1) and (int(index1.split("_")[-1]), int(index2.split("_")[-1])) in self.index_list:
                    if (index1, index2) not in index_blast_dict:
                        index_blast_dict[(index1, index2)] = {} 
                    else:
                        pass 
                    #print(index1,index2)
                    index_blast_dict[(index1, index2)][target] = {} 
                    index_blast_dict[(index1, index2)][target]["R1"]   = None
                    index_blast_dict[(index1, index2)][target]["R2"]   = None
                    index_blast_dict[(index1, index2)][target]["R1"]   = R1_hit_dict 
                    index_blast_dict[(index1, index2)][target]["R2"]   = R2_hit_dict
            #print(index_blast_dict.keys())
            return index_blast_dict
        
                 
        def single_process(fasta,_id):
            os.chdir("tmp")
            os.mkdir(str(_id))
            os.chdir(str(_id)) 
            
            fastaR1 = fasta 
            fastaR2 = fasta.replace("R1","R2") 
            blastR1 = fastaR1 + ".out" 
            blastR2 = fastaR2 + ".out"
            name = fastaR1.replace(".fasta","").replace("R1_","").split("/")[-1]
            
            R1_fasta_dict    = read_fasta("../../fragmented_fasta/" + fastaR1)
            R1_fastq_dict    = read_fastq("../../fragmented_fastq/" + fastaR1.replace("fasta","fastq")) 
            R2_fasta_dict    = read_fasta("../../fragmented_fasta/" + fastaR2)
            R2_fastq_dict    = read_fastq("../../fragmented_fastq/" + fastaR2.replace("fasta","fastq")) 
            fasta_dict = {"R1":{}, "R2":{}}
            fastq_dict = {"R1":{}, "R2":{}}
            blast_dict = {} 
            for key in R1_fasta_dict:
                blast_dict[key] = {} 
                blast_dict[key]["R1"] = {}  
                blast_dict[key]["R2"] = {}
                blast_dict[key]["R1"]["seq"] = R1_fasta_dict[key] 
                blast_dict[key]["R2"]["seq"] = R2_fasta_dict[key]
                blast_dict[key]["R1"]["hit"] = collections.defaultdict(list)
                blast_dict[key]["R2"]["hit"] = collections.defaultdict(list)
            blast_dict       = read_blast("../../blast_result/" + fastaR1 + ".out", blast_dict, direction="R1") 
            blast_dict       = read_blast("../../blast_result/" + fastaR2 + ".out", blast_dict, direction="R2")
            index_blast_dict = search_index(self.index_dict, blast_dict)
          
            for index_set in index_blast_dict:
                index_dir = "{}_{}".format(index_set[0],index_set[1])
                os.mkdir(index_dir)
                with open("{}/{}_R1.fasta".format(index_dir,index_dir),"w") as fasta:
                    for read in index_blast_dict[index_set]:
                        fasta.write(">{}\n".format(read))
                        fasta.write(R1_fasta_dict[read] + "\n") 
                
                with open("{}/{}_R2.fasta".format(index_dir,index_dir),"w") as fasta:
                    for read in index_blast_dict[index_set]:
                        fasta.write(">{}\n".format(read))
                        fasta.write(R2_fasta_dict[read] + "\n") 
                
                index_nums = (int(index_set[0].split("_")[-1]), int(index_set[1].split("_")[-1]))
                
                flag1, flag2 = 0, 0
                
                fasta_file_dict = {} 
                fastq_file_dict = {} 
                if index_nums[0] != 0:
                    keys_array_R1 = []
                    with open("./{}/FW.fasta".format(index_dir),"w") as o:
                        for project_id, sample_id, primer in self.index_primer_dict[index_nums]["FW"]:
                            keys = [] 
                            o.write(">FW;{};{}\n".format(project_id,sample_id)) 
                            o.write(primer + "\n") 
                            fasta_file_dict[(project_id,sample_id)] = open("{}_{}_R1.fasta".format(project_id,sample_id), "w")
                            fastq_file_dict[(project_id,sample_id)] = open("{}_{}_R1.fastq".format(project_id,sample_id), "w")

                    os.system("blastn -query {}/{}_R1.fasta -subject {}/FW.fasta -task blastn-short -evalue 1.0e-3 -outfmt \"6 sseqid qseqid\" -out {}_FW.out".format(index_dir,index_dir,index_dir,index_dir))
                    with open(str(index_dir) + "_FW.out") as f:
                        lines = f.readlines() 
                        for line in lines:
                            sseqid, key = line.rstrip().split("\t") 
                            project_id  = sseqid.split(";")[1] 
                            sample_id   = sseqid.split(";")[2] 
                            fasta = fasta_file_dict[(project_id,sample_id)]
                            fasta.write(">{}\n".format(key))
                            fasta.write(R1_fasta_dict[key] + "\n") 
                            
                            fastq = fastq_file_dict[(project_id,sample_id)]
                            fastq.write("@{}\n".format(R1_fastq_dict[key][0]))
                            fastq.write(R1_fastq_dict[key][1] + "\n") 
                            fastq.write(R1_fastq_dict[key][2] + "\n") 
                            fastq.write(R1_fastq_dict[key][3] + "\n") 
                    for key in fastq_file_dict:
                        fastq_file_dict[key].close()
                    for key in fasta_file_dict:
                        fasta_file_dict[key].close()
                    flag1 = 1
                
                fasta_file_dict = {} 
                fastq_file_dict = {} 
                if index_nums[1] != 0:
                    with open("./{}/RV.fasta".format(index_dir),"w") as o:
                        for project_id, sample_id, primer in self.index_primer_dict[index_nums]["RV"]:
                            keys = []
                            o.write(">RV;{};{}\n".format(project_id,sample_id))
                            o.write(primer + "\n") 
                            fasta_file_dict[(project_id,sample_id)] = open("{}_{}_R2.fasta".format(project_id,sample_id), "w")
                            fastq_file_dict[(project_id,sample_id)] = open("{}_{}_R2.fastq".format(project_id,sample_id), "w")

                    os.system("blastn -query {}/{}_R2.fasta -subject {}/RV.fasta -task blastn-short -evalue 1.0e-3 -outfmt \"6 sseqid qseqid\" -out {}_RV.out".format(index_dir,index_dir,index_dir,index_dir))
                    with open(str(index_dir) + "_RV.out") as f:
                        lines = f.readlines()
                        for line in lines:
                            sseqid, key = line.rstrip().split("\t") 
                            project_id  = sseqid.split(";")[1] 
                            sample_id   = sseqid.split(";")[2] 
                            fasta = fasta_file_dict[(project_id,sample_id)]
                            fasta.write(">{}\n".format(key))
                            fasta.write(R2_fasta_dict[key] + "\n") 
                            
                            fastq = fastq_file_dict[(project_id,sample_id)]
                            fastq.write("@{}\n".format(R1_fastq_dict[key][0]))
                            fastq.write(R2_fastq_dict[key][1] + "\n") 
                            fastq.write(R2_fastq_dict[key][2] + "\n") 
                            fastq.write(R2_fastq_dict[key][3] + "\n") 

                    for key in fastq_file_dict:
                        fastq_file_dict[key].close()
                    for key in fasta_file_dict:
                        fasta_file_dict[key].close()
                    flag2 = 1 
                
            os.chdir("../../") 
        
        if "blast_result" not in os.listdir("."):
            os.mkdir("./blast_result")
        
        if "demultiplexed" not in os.listdir("."): 
            os.mkdir("demultiplexed") 
        
        os.chdir("demultiplexed")
        for project_id, sample_id in self.experimental_dict:
            if project_id not in os.listdir("."):
                os.mkdir(str(project_id))
            os.chdir(str(project_id))
            if sample_id not in os.listdir("."):
                os.mkdir(str(sample_id))
            os.chdir(str(sample_id)) 
            os.system("touch  R1.fasta R2.fasta R1.fastq R2.fastq")
            os.chdir("../../")
        os.chdir(pwd) 
        
        os.chdir("./fragmented_fasta")
        if sge:
            fastas = [f for f in os.listdir(".") if ".fasta" in f and ".out" not in f]
            fastas = [fastas[i:i+50] for i in range(0,len(fastas),50)] 
            for f in fastas:
                f = str(tuple(f))
                f = f.replace("'","").replace("(","").replace(")","").replace(",","")
                os.system("qsub" + " -o ./output/ " + " -e ./error/ " + " -v " + "TARGET=\"" + f + "\" ./script/blast.sh")
        else:
            commands = []
            for fasta in self.fasta_list:        
                #print(fasta)
                commands.append("blastn -query {} -db ../primers/primers -task blastn-short -evalue 1.0e-3 -outfmt \"6 qseqid qseq sseqid sseq sstrand pident qlen length mismatch gapopen qstart qend sstart send gaps evalue bitscore btop\" -out ../blast_result/{}.out".format(fasta,fasta))    
                if len(commands) % num_threads == 0:
                    os.system("&".join(commands)) 
                    commands = []   
                else:
                    pass 

            if len(commands) > 0:
                os.system("&".join(commands)) 
        os.chdir(pwd)
        
        #os.chdir("demultiplexed") 
        if "tmp" not in os.listdir("."):
            os.mkdir("tmp")
 
        ps    = []
        for process_id, fasta in enumerate([f for f in os.listdir("./fragmented_fasta") if "R1" in f]):
            ps.append(mp.Process(target=single_process,args=(fasta,process_id)))
        
        finishedList = []               
        presentindex = 0
        for p in ps[0:24]:
            p.start()
            presentindex += 1

        #print "start"
        collected = 0
        while 1:
            time.sleep(1)
            for i, p in enumerate(ps[:presentindex]):
                if p.is_alive() or i in finishedList:
                    if p.is_alive():
                        collected += 1
                else:
                    print(i) 
                    finishedList.append(i)                   
                    if presentindex < len(ps):
                        ps[presentindex].start()
                        presentindex += 1
            
            if len(finishedList) == len(ps):
                break
        
        for i in os.listdir("tmp/"):
            commands = []
            os.chdir("tmp/{}".format(i)) 
            for project_id, sample_id in self.experimental_dict:
                if "{}_{}_R1.fasta".format(project_id,sample_id) in os.listdir("./"):
                    commands.append("cat {}_{}_R1.fasta >> ../../demultiplexed/{}/{}/R1.fasta".format(project_id,sample_id,project_id,sample_id)) 
                if "{}_{}_R2.fasta".format(project_id,sample_id) in os.listdir("./"):
                    commands.append("cat {}_{}_R2.fasta >> ../../demultiplexed/{}/{}/R2.fasta".format(project_id,sample_id,project_id,sample_id)) 
                if "{}_{}_R1.fastq".format(project_id,sample_id) in os.listdir("./"):
                    commands.append("cat {}_{}_R1.fastq >> ../../demultiplexed/{}/{}/R1.fastq".format(project_id,sample_id,project_id,sample_id)) 
                if "{}_{}_R2.fastq".format(project_id,sample_id) in os.listdir("./"):
                    commands.append("cat {}_{}_R2.fastq >> ../../demultiplexed/{}/{}/R2.fastq".format(project_id,sample_id,project_id,sample_id)) 
            for j in range(0,len(commands),32):
                os.system(" & ".join(commands[j:j+32]))
            os.chdir("../../") 
            os.system("rm -rf tmp/{}".format(i))
        os.chdir(pwd)
    
    def save_sample_state(self):
        for project_id, sample_id in self.experimental_dict:
            with open("demultiplexed/{}/{}/sample.json".format(project_id,sample_id),"w") as o:
                json.dump(self.experimental_dict[(project_id,sample_id)].__dict__,o)
                
 
def process(anaysis,project_id,sample_id,target,min_match=0.8,read_type="merged"):
    os.chdir("./" + project_id + "/" + sample_id)
    sample = analysis.experimental_dict[(project_id,sample_id)]
    sample.min_match = min_match
    if sample.control[1] == "control" or sample.control[1] == "N.A.":
        pass 
    else:
        sample.conobj = analysis.experimental_dict[sample.control]
    sample.merge()
    #sample.alignment(read_type,num_threads=1)
    #sample.extract_target_regions(read_type)
    #if sample.method == "barcode_count":
    #    sample.barcode_count(read_type)
    #else:
    #    sample.calc_editing_spectrum(read_type, be=True) 
    #sample.alignment("R1",num_threads=1)
    #sample.extract_target_regions("R1")
    #sample.calc_editing_spectrum("R1", be=True) 
    #sample.alignment("R2",num_threads=1)
    #sample.extract_target_regions("R2")
    #sample.calc_editing_spectrum("R2", be=True) 

    os.chdir("../../")
    
def parallel(ps): 
    finishedList = []               
    presentindex = 0
    for p in ps[0:mp.cpu_count()]:
        p.start()
        presentindex += 1
    collected = 0
    while 1:
        time.sleep(1)
        for i, p in enumerate(ps[:presentindex]):
            if p.is_alive() or i in finishedList:
                if p.is_alive():
                    collected += 1
            else:
                print(i)
                finishedList.append(i)
                if presentindex < len(ps):
                    ps[presentindex].start()
                    presentindex += 1
        
        if len(finishedList) == len(ps):
            break

if __name__ == "__main__":
    pwd = os.getcwd()
    print(pwd) 
    p = argparse.ArgumentParser()
    p.add_argument("-s","--sample_csv",  type=str,   default   = "None",   help=" ")
    p.add_argument("-c","--const_fasta", type=str,   default   = "None",   help=" ")
    p.add_argument("-m","--miseq_data",  type=str,   default   = "None",   help=" ")
    p.add_argument("-i","--min_match",   type=float, default   = 0.5,      help=" ")
    p.add_argument("-t","--target",      type=str,   default   = "merged", help=" ")
    
    args   = p.parse_args()
    target = args.target
    analysis = ANALYSIS() 
    analysis.load_setting(args.sample_csv,args.const_fasta) 
    print("Loading fastq and spliting") 
    analysis.load_data(args.miseq_data) 
    print("Converting fastq to fasta") 
    analysis.fastq2fasta() 
    print("Searching index and primer sequences") 
    os.chdir(pwd)
    analysis.search_primer() 
    os.chdir(pwd)
    analysis.save_sample_state() 
    os.chdir("demultiplexed")     
    ps = []
    process_id = 0
    """
    for project_id, sample_id in analysis.experimental_dict:
        if sample_id not in os.listdir("./" + project_id + "/"):
            pass
        elif analysis.experimental_dict[(project_id,sample_id)]["control"][1] == "control":
            ps.append(mp.Process(target=process,args=(analysis,project_id,sample_id,process_id,args.min_match,args.target)))
            process_id += 1
    print("Visualization of control samples")
    parallel(ps)
        
    ps = []
    process_id = 0
    for project_id, sample_id in analysis.experimental_dict:
        if sample_id not in os.listdir("./" + project_id + "/"):
            pass
        elif analysis.experimental_dict[(project_id,sample_id)]["control"] != "control":
            ps.append(mp.Process(target=process,args=(analysis,project_id,sample_id,process_id,args.min_match,args.target)))
            process_id += 1 
    parallel(ps)
    """
