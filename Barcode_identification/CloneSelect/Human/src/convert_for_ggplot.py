import collections
from commandr import command, Run


@command('matrix_to_long_format')
def matrix_to_long_format(mat_file):
    bc_name_list = []
    with open(mat_file, 'r') as fin:
        for line in fin:
            if line.startswith('clone_name'):
                c = line.strip('\n').split('\t')
                bc_name_list = c

    bc_name_list = bc_name_list[1:]  # remove matrix header

    is_match = 'NT'

    print('order,query_bc,norm_bc,bc,is_match')
    order = 1
    with open(mat_file, 'r') as fin:
        for line in fin:
            if line.startswith('clone_name'):
                continue
            col = line.strip('\n').split('\t')
            query_bc_name = col[0]

            for bc_value, bc_name in zip(col[1:], bc_name_list):
                if query_bc_name == bc_name:
                    is_match = 'OT'
                print('%d,%s,%s,%s,%s' %
                      (order, query_bc_name, bc_value, bc_name, is_match))
                order += 1
                is_match = 'NT'


if __name__ == '__main__':
    Run()
