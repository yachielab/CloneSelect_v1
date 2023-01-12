import sys
import collections
import random


if __name__ == '__main__':
    # all barcode types (main figure)
    bc_matrix = 'rep1/csv/barcode_counts_rep1_normalized.csv'
    # only AUG converted barcodes
    bc_matrix = 'rep1/csv/barcode_counts_rep1_only_atg_normalized.csv'

    bc_name_list = []
    with open(bc_matrix, 'r') as fin:
        for line in fin:
            if line.startswith('Barcode'):
                c = line.strip('\n').split(',')
                bc_name_list = c

    bc_name_list = bc_name_list[1:]  # remove matrix header

    barcode_name_map = {}
    data = collections.defaultdict(str)
    is_match = 'NT'

    print('order,query_bc,norm_bc,bc,is_match')
    order = 1
    with open(bc_matrix, 'r') as fin:
        for line in fin:
            if line.startswith('Barcode'):
                continue
            col = line.strip('\n').split(',')
            query_bc_name = col[0]

            for bc_value, bc_name in zip(col[1:], bc_name_list):
                if query_bc_name == bc_name:
                    is_match = 'OT'
                print('%d,%s,%s,%s,%s' %
                      (order, query_bc_name, bc_value, bc_name, is_match))
                order += 1
                is_match = 'NT'
