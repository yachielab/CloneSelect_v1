import os 
import re
import sys 
import copy
sys.path.append("/home/t12870hm/library/CLOVER/")
sys.path.append("/home/t12870hm/library/CLOVER/CLOVER")
from CLOVER import analyzer
from CLOVER import integrator
from CLOVER.visualizer import *
import matplotlib.pyplot as plt

targets  = ["112-V3-BC8","112-V3-BC9","112-V3-BC16","112-V3-BC19","112-V9-BC18","112-V9-BC19","112-V9-BC21","112-V10-BC14","112-V10-BC15","112-V10-BC17","112-V10-BC25","112-V10-BC29","112-V10-BC31","112-V10-BC33","112-V4-BC7","112-V4-BC8"] 
template = {
            "project":       None,
            "reference_seq": "CTACTAGAGGATCTATTTCCGGTGAATTCCCGAGCGTGTCAGGGTGACCGTGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGG",
            "guide_seq":     "TCCGGTGAATTCCCG",
            "zero_position": ("TCCGGTGAATTCCCG", 15),
            "call_range":    (-25, 50),
            "read1":         "/home/t12870hm/MiSeq/MiSeq_package/cloneselect2022/cloneselect_barcode_count/miseq_data/R1/{}_Rep1_R1.fastq",
            "read2":         "/home/t12870hm/MiSeq/MiSeq_package/cloneselect2022/cloneselect_barcode_count/miseq_data/R2/{}_Rep1_R2.fastq",
            "attribute":     None,
            "output_path":   "/home/t12870hm/MiSeq/MiSeq_package/cloneselect2022/cloneselect_barcode_count/output/Rep1/{}",
            "min_identity":  0.8
           }

samples  = [] 
Analyzer = analyzer.Analyzer
for target in targets:
    sample_info                = copy.deepcopy(template) 
    sample_info["project"]     = target
    sample_info["attribute"]   = {"project": target}
    sample_info["read1"]       = sample_info["read1"].format(target)
    sample_info["read2"]       = sample_info["read2"].format(target)
    sample_info["output_path"] = sample_info["output_path"].format(target)
    os.makedirs(sample_info["output_path"], exist_ok=True)
    sample = Analyzer(**sample_info)
    samples.append(sample) 

Integrator  = integrator.Integrator
int_samples = Integrator(*samples)
#int_samples.pre_processing(pickle=True, cpunum=4)
int_samples.load_pickle()
barcodes_dict = int_samples.get_specificpattern(query="TTCCGGTGAATTCCC[ATGC]{14,20}[ATGC]TGGTGAGCAAGGGC", call_range=(-25, 50), style="sample", cpunum=16)

for project in targets:
    for key in list(barcodes_dict[project]["pattern"].keys()):
        match = re.search("TTCCGGTGAATTCCC([ATGC]{14,20}[ATGC]TG)GTGAGCAAGGGC", key)
        print(project, key, match.group(1), barcodes_dict[project]["pattern"][key]["count"], sep="\t")

