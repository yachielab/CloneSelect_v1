#!/usr/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N reamp_lib
#read1 contains cellBC, UMI; read2 contains transcript reads

fwfastq=$1
rvfastq=$2


#trim 5' end index and adaptor for read2
cutadapt \
  -e 0.1 --no-indels -O 32 \
  -G $EGFP \
  --pair-filter=both \
  --discard-untrimmed \
  -o out.1.fastq -p out.2.fastq \
  $fwfastq $rvfastq

#Trim 3' end
cutadapt \
  -e 0.1 --no-indels -O 20\
  -A TCTAGAGTATTACCATGGTC \
  --discard-untrimmed \
  --pair-filter=both \
  -o reads.1_1.fastq -p reads.2_1.fastq \
  out.1.fastq out.2.fastq

#merge the dnBC(read2) and cellBC (read1)
awk '{ getline other <"reads.1_1.fastq" } { print $0 (FNR % 2 == 0 ? other : "") }' reads.2_1.fastq > BClibmerge.fastq
wait

echo ' === reads merged as : "BClibmerge.fastq"'

#bartender barcode extraction
##output: BClibmerge_extracted
##q-score above 30
##length of BC 50. G[50]T
##no mismatch allowed
bartender_extractor_com \
     -f BClibmerge.fastq \
     -o BClibmerge_extracted \
     -p  AGTAA[50]TTTTT \
     -q ? \
     -m 0
wait

#split the merged BC
awk -F, 'FNR > 1 {print substr($1,1,30) "," substr($1,31,50) "," $2}' BClibmerge_extracted_barcode.txt > BClibsplit.csv
