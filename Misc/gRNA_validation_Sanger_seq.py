import os

def name2well(mapping_file):
    data = {}
    with open(mapping_file, 'r') as fh:
        for line in fh:
            cols = line.strip().split(',')
            well, sample = cols[0], cols[1]
            data[sample] = well
    return data

def load_BC_reference(bc_file):
    data = {}
    with open(bc_file, 'r') as fh:
        for line in fh:
            cols = line.strip().split(',')
            bc_id, bc_seq, well_id = cols[0], cols[1], cols[3]
            data[well_id] = {'bc_seq': bc_seq, 'bc_id': bc_id}
    return data

    
def reverse_comp(x):
    return ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
    

if __name__ == '__main__':
    well_map = name2well('well_sanger_map.csv')
    reference_bc_data = load_BC_reference('72_bcid2bcseq_ref.csv')
    
    seq_list = filter(lambda x: x.endswith('.txt'), os.listdir('sanger_data'))
    data = {}
    for seq_file in seq_list:
        sample_name = seq_file.split('.')[0]
        seq = ''
        with open(os.path.join('sanger_data', seq_file), 'r') as fh:
            for line in fh:
                seq += line.strip()

        well_id = well_map[sample_name]
        a, b = well_id[0], well_id[1:]
        data[sample_name] = {'well_id': well_id, 'sample_name': sample_name,
                             'seq': seq, 'a': a, 'b': b}
        #print sample_name, well_map[sample_name]

    # sorted by well id
    scaffold = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTA'
    for k, v in sorted(data.items(), key=lambda x: (x[1]['a'], x[1]['b'])):
        pos = v['seq'].find(scaffold) # found scaffold with perfect match
        if pos != -1:
            start = pos-20
            end = pos
            sgRNA_seq = v['seq'][start:end]
            sample_well_id = v['well_id']
            candidate = reverse_comp(sgRNA_seq)
            ref = reference_bc_data[sample_well_id]['bc_seq']
            if candidate == ref:
                #print k, v['well_id'], 'Confirmed'
                print '>%s\n%s' % (k, v['seq'])
            else:
                pass #print k, v['well_id'], 'sgRNA_not_matched_ref_BC'
            
        elif pos == -1:
            pass#print k, v['well_id'], 'error'
            
            



        
    
