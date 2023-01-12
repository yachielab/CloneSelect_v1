import collections
from commandr import command, Run

MIN_BARCODE_LEN = 16


@command('normalize_barcodes')
def normalize_barcodes(starcode_input):
    # input should be starcode output

    barcode_db = collections.defaultdict(str)
    with open(starcode_input, 'r') as fin:
        for line in fin:
            cols = line.strip().split('\t')
            barcode = cols[0]

            if len(barcode) >= MIN_BARCODE_LEN:
                count = int(cols[1])
                barcode_db[barcode] = count

    total_count = sum(barcode_db.values())
    rank = 1

    print('barcode,norm_freq,RPM,raw_count')
    for k, v in sorted(barcode_db.items(), key=lambda x: x[1], reverse=True):
        norm_freq = v/total_count
        rpm = v*1e6/total_count
        print("%s,%f,%f,%d" % (k, norm_freq, rpm, v))

        rank += 1


if __name__ == '__main__':
    Run()
