import os 
path1   = "/home/soh.i/illumina/02272018_ClonSelect_MiSeq/data/production/refined/R1/"
path2   = "/home/soh.i/illumina/02272018_ClonSelect_MiSeq/data/production/refined/R2/"
targets = ["112-V3-BC8","112-V3-BC9","112-V3-BC16","112-V3-BC19","112-V9-BC18","112-V9-BC19","112-V9-BC21","112-V10-BC14","112-V10-BC15","112-V10-BC17","112-V10-BC25","112-V10-BC29","112-V10-BC31","112-V10-BC33","112-V4-BC7","112-V4-BC8"] 
suffix1 = "_Rep2_R1.fastq"
suffix2 = "_Rep2_R2.fastq"
targets.sort()
for f in targets:
    f1 = f + suffix1
    f2 = f + suffix2
    #print "python count_barcode.py {} {}".format(path + f,"112-96lib-identified-known-barcode-gDNA.csv")
    #os.system("python count_barcode.py {} {}".format(path + f,"112-96lib-identified-known-barcode-gDNA.csv")) 
    print "python2 ./script/count_barcode.py {} {} {}".format(path1 + f1, path2 + f2, "112-96lib-identified-known-barcode-gDNA.csv")
    os.system("python2 ./script/count_barcode.py {} {} {}".format(path1 + f1, path2 + f2, "112-96lib-identified-known-barcode-gDNA.csv")) 

