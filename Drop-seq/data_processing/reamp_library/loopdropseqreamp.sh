#!/usr/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N demultiplex
#to run dropseqreamp.sh for mutiple samples

source ~/.bashrc
bash

for f in ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/061622_Hiseq_dnBCReamp/demultiplex/*.1.fastq.gz;
do
    cd ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/061622_Hiseq_dnBCReamp/library
    filename=${f##*/}
    filename=${filename%%.*}
    mkdir -p ./${filename}
    cd ./${filename}
    wait
    bash ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/061622_Hiseq_dnBCReamp/library/dropseqreamp.sh \
    ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/061622_Hiseq_dnBCReamp/demultiplex/${filename}.1.fastq.gz \
    ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/061622_Hiseq_dnBCReamp/demultiplex/${filename}.2.fastq.gz
    echo "running dropseq reamp on ${filename}.1.fastq.gz ${filename}.2.fastq.gz"
done


#output number of barcodes per index
for f in ../demultiplex_Hiseq_HKTNNBCX3/*.1.fastq.gz;
do
    cd ~/DryLab/Projects/scCloneSelect/scCSExp17.3-3/061622_Hiseq_dnBCReamp/library
    filename=${f##*/}
    filename=${filename%%.*}
    cd ./${filename}
    wait
    count=`echo $(cat "BClibsplit.csv"|wc -l)`
    wait
    echo "$filename has $count barcodes"
done > counts.txt
