import collections
import csv
from commandr import command, Run


def generate_db(_file):
    db = collections.defaultdict(str)
    with open(_file, 'r') as fin:
        reader = csv.DictReader(fin)
        for cols in reader:
            db[cols['barcode']] = cols
    return(db)


@command('merge_pDNA_gDNA')
def merge_pDNA_gDNA(pDNA_file, gDNA_file):
    db_pDNA = generate_db(pDNA_file)
    db_gDNA = generate_db(gDNA_file)

    # get unique set of barcodes from two datasets
    bcs = set()
    for v in db_pDNA.values():
        bcs.add(v['barcode'])
    for v in db_gDNA.values():
        bcs.add(v['barcode'])

    new_db = collections.defaultdict(str)
    for barcode in bcs:
        pdna = db_pDNA.get(barcode)
        gdna = db_gDNA.get(barcode)
        if pdna and gdna:
            new_db[barcode] = {'norm_pfreq': pdna["norm_freq"], 'pRPM': pdna["RPM"], 'p_raw': pdna["raw_count"],
                               'norm_gfreq': gdna["norm_freq"], 'gRPM': gdna["RPM"], 'g_raw': gdna["raw_count"]}

    _id = 1

    print('bc,clone_name,pfreq,pRPM,p_raw,gfreq,gRPM,g_raw')
    for barcode, v in sorted(new_db.items(), key=lambda x: x[1]['norm_gfreq'], reverse=True):
        clone_name = 'Clone %d' % (_id)
        _id += 1
        print("%s,%s,%s,%s,%s,%s,%s,%s" % (
            barcode+"GTG", clone_name, v["norm_pfreq"], v["pRPM"], v["p_raw"], v["norm_gfreq"], v["gRPM"], v["g_raw"]))


if __name__ == '__main__':
    Run()
