#!/bin/sh

starcode="/Users/leon/dev/starcode/starcode"
${starcode} -i data/112_96_pDNA_barcodes_rev_conv.tsv --dist 4 -o data/112_96_pDNA_barcodes_starcode.csv  --print-clusters
${starcode} -i data/112_96_gDNA_barcodes_rev_conv.tsv --dist 4 -o data/112_96_gDNA_barcodes_starcode.csv  --print-clusters
