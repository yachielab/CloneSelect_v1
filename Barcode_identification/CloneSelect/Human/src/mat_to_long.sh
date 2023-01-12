#!/bin/sh

# Wrapper for src/convert_for_ggplot.py

python="/Users/leon/.pyenv/shims/python"

# Rep1 
${python} src/convert_for_ggplot.py matrix_to_long_format -m data/GFP_isolation/rep1_result/barcode_counts_rep1_normalized.tsv > barcode_counts_rep1_normalized_long.csv 
${python} src/convert_for_ggplot.py matrix_to_long_format -m data/GFP_isolation/rep1_result/barcode_counts_rep1_only_atg_normalized.tsv > barcode_counts_rep1_atg_normalized_long.csv

# Rep2
${python} src/convert_for_ggplot.py matrix_to_long_format -m data/GFP_isolation/rep2_result/barcode_counts_rep2_normalized.tsv > barcode_counts_rep2_normalized_long.csv 
${python} src/convert_for_ggplot.py matrix_to_long_format -m data/GFP_isolation/rep2_result/barcode_counts_rep2_only_atg_normalized.tsv > barcode_counts_rep2_atg_normalized_long.csv

