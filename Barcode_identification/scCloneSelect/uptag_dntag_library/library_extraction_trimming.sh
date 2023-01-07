#!/bin/bash

#scCSlibrary3.sh
#run by bash barcode_extraction_trimming.sh Read1_fastq Read2_fastq

source ~/.bashrc

fwfastq=$1
rvfastq=$2


# Filter and trim reads by defining the sequence 20 bp before and after barcodes
cutadapt \
    -e 0.1 --no-indels -O 20\
    -o BClib.1.fastq.gz -p BClib.2.fastq.gz \
    -g TCTATTTCCGGTGAATTCCC...GTGGTGAGCAAGGGCGAGGA \
    -G GACCATGGTAATACTCTAGA...TACTTGTACAGCTCGTCCAT \
    --discard-untrimmed \
    $fwfastq $rvfastq
wait

echo '### reads extracted and trimed'

#unzip files
gunzip BClib.1.fastq.gz
gunzip BClib.2.fastq.gz

#merge the fw and reverse reads
awk '{ getline other <"BClib.2.fastq" } { print $0 (FNR % 2 == 0 ? other : "") }' BClib.1.fastq > BClibmerge.fastq
wait

echo ' === reads merged as : "BClibmerge.fastq"'

#bartender barcode extraction
##output: BClibmerge_extracted
##q-score above 30
##length of BC 47. G[47]T
##no mismatch allowed
bartender_extractor_com \
     -f BClibmerge.fastq \
     -o BClibmerge_extracted \
     -q ? \
     -p  G[47]T \
     -m 0
wait

#bartender clustering
bartender_single_com \
    -f BClibmerge_extracted_barcode.txt \
    -o BClibmerge \
    -d 5 \
    -z -1 \
    -l 3 \
    -s 1
wait
echo '### bartender clustering complete"'

#split the merged BC
awk -F, 'FNR > 1 {print substr($2,1,30) "," substr($2,31,50) "," $4}' BClibmerge_cluster.csv > BClibsplit.csv
wait
echo '### split bartender output saved as "BClibsplit.csv"'

#proceeds to running the barcode filtering and quantification by python
python ../library_filter_quantification.py BClibsplit.csv
