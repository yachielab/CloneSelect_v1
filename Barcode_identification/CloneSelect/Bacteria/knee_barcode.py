import collections
from os import sep
from kneed import DataGenerator, KneeLocator
from collections import defaultdict
from commandr import command, Run
import csv


def get_knee_rank(count_mat):
    rank = list(range(0, len(count_mat)))
    i = 1
    cum_count = []
    for _, v in sorted(count_mat.items(), key=lambda x: x[1], reverse=True):
        bc_count = sum(v)
        i += bc_count
        cum_count.append(i)

    kneedle = KneeLocator(
        rank, cum_count, S=1, curve='concave', direction='increasing')
    return kneedle.knee


@command('filter_bc_by_knee')
def filter_bc_by_knee(count_file):
    '''
    count_file: raw count barcode csv.
    '''

    count_mat = defaultdict(str)
    tot_rep1 = 0
    tot_rep2 = 0
    tot_rep3 = 0
    with open(count_file) as fin:
        for line in fin:
            cols = line.strip('\n').split('\t')
            if 'barcode' in cols[0]:
                continue
            barcode = cols[0]
            rep1, rep2, rep3 = int(cols[1]), int(cols[2]), int(cols[3])
            count_mat[barcode] = {'rep1': rep1, 'rep2': rep2,
                                  'rep3': rep3}
            tot_rep1 += rep1
            tot_rep2 += rep2
            tot_rep3 += rep3

    min_barcode_count = 0

    barcode_filt = collections.defaultdict(str)
    mean_barcode_count = collections.defaultdict(list)
    for k, d in sorted(count_mat.items()):
        if d['rep1'] > min_barcode_count and d['rep2'] > min_barcode_count and d['rep3'] > min_barcode_count:
            mean_raw_count = sum([d['rep1'], d['rep2'], d['rep3']])/3.0

            # Barcode freq
            freq_rep1 = d['rep1']/tot_rep1
            freq_rep2 = d['rep2']/tot_rep2
            freq_rep3 = d['rep3']/tot_rep3
            freq_mean = sum([freq_rep1, freq_rep2, freq_rep3])/3.0

            # RPM barcode
            RPM_rep1 = d['rep1']*1e6/tot_rep1
            RPM_rep2 = d['rep2']*1e6/tot_rep2
            RPM_rep3 = d['rep3']*1e6/tot_rep3
            RPM_mean = sum([RPM_rep1, RPM_rep2, RPM_rep3])/3.0

            # print(k, d, mean_raw_count, freq_mean, RPM_mean)
            barcode_filt[k] = {
                'rep1': d['rep1'], 'rep2': d['rep2'], 'rep3': d['rep3'],
                'rep1_RPM': RPM_rep1, 'rep2_RPM': RPM_rep2, 'rep3_RPM': RPM_rep3,
                'rep1_freq': freq_rep1, 'rep2_freq': freq_rep2, 'rep3_freq': freq_rep3,
                'mean_raw_count': mean_raw_count,
                'mean_freq': freq_mean, 'RPM_mean': RPM_mean}
            mean_barcode_count[k].append(RPM_mean)

    knee = get_knee_rank(mean_barcode_count)
    rank = 1

    print('rank,knee_flag,barcode,rep1,rep2,rep3,RPM_r1,RPM_r2,RPM_r3,freq_r1,freq_r2,freq_r3,mean_raw_count,mean_freq,RPM_mean')
    for k, v in sorted(barcode_filt.items(), key=lambda x: x[1]['RPM_mean'], reverse=True):
        if rank <= knee:
            knee_flag = 'yes'
        else:
            knee_flag = 'no'

        print(
            f"{rank},{knee_flag},\"{k}\",{v['rep1']}, {v['rep2']}, {v['rep3']}, {v['rep1_RPM']}, {v['rep2_RPM']}, {v['rep3_RPM']}, {v['rep1_freq']}, {v['rep2_freq']}, {v['rep3_freq']}, {v['mean_raw_count']}, {v['mean_freq']}, {v['RPM_mean']}")
        rank += 1


if __name__ == '__main__':
    Run()
