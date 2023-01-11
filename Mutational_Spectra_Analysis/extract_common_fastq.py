import os 
import sys 

def read_fastq(fastq_name):
    seq_dict = {} 
    with open(fastq_name) as f:
        n = 0 
        for line in f:
            if line[0] == "@" and n % 4 == 0:
                key = line[1:].rstrip() 
                key = key.split(" ")[0] 
                key = key.replace(":","_")
                seq_dict[key] = [line.rstrip()] 
            else:
                seq_dict[key].append(line.rstrip()) 
            n += 1
    return seq_dict 

def merge(f1, f2):   
    R1_dict = read_fastq(f1) 
    R2_dict = read_fastq(f2) 
    with open("R1_common.fastq","w") as R1_tmp, open("R2_common.fastq","w") as R2_tmp:
        for key in list(set(R1_dict.keys()) & set(R2_dict.keys())):
            for i in range(4):
                R1_tmp.write(R1_dict[key][i] + "\n") 
                R2_tmp.write(R2_dict[key][i] + "\n") 

if __name__ == "__main__":
    extract(sys.argv[1], sys.argv[2]) 

