import csv
import collections

if __name__ == '__main__':
    # load zeocin colony count data
    plate_df = collections.defaultdict(int)
    plate_data = 'data/ECS_BC100pool_heatmap_data_08052021.tsv'
    with open(plate_data, 'r') as fin:
        reader = csv.DictReader(fin, delimiter='\t')
        for row in reader:
            plate_df[row['id']] = row

    # load barcode count data
    barcode_count = "../Barcode_rank_plot/data/BC100_processed_05102022.csv"
    barcode_df = collections.defaultdict(str)
    with open(barcode_count, 'r') as fin:
        reader = csv.DictReader(fin, delimiter=',')
        for row in reader:
            barcode_df[row['barcode']] = row

    for _, plate_val in plate_df.items():
        bcseq = plate_val['BC seq detected']
        if bcseq[:-3] in barcode_df.keys():
            k = bcseq[:-3]
            bc_freq = barcode_df[k]['ave.freq']
            print(plate_val['id'], plate_val['no'], plate_val['BC'], plate_val['plate'], plate_val['plate'], plate_val['BC seq detected'],
                  plate_val['is_expected_BC'], plate_val['is_CAA_conv'], bc_freq, plate_val['Total_colony'], plate_val['Total_Zeo_colony'], plate_val['Zeo_colony_percent'])
