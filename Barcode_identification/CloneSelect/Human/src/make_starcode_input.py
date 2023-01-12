import commandr


@commandr.command('make_starcode')
def make_starcode(count_data):
    with open(count_data) as fin:
        for line in fin:
            if not line.startswith('Reference'):
                cols = line.strip().split('\t')
                barcode = cols[2]
                count = cols[3]
                print("%s\t%s" % (barcode, count))


if __name__ == '__main__':
    commandr.Run()
