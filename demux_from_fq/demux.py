#!/usr/bin/python

__author__  = 'Shan Sabri'
__email__   = 'ShanASabri@gmail.com'
__date__    = '2019/09/11'
__version__ = '1.1'


import datetime, gzip

def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}


def write_record(r, out):
    out.write("%s\n%s\n%s\n%s\n" % (r["name"], r["sequence"], r["optional"], r["quality"]))

def demux(fq, barcodes):
    o1 = gzip.open(barcodes[0] + ".fastq.gz", 'w')
    o2 = gzip.open(barcodes[1] + ".fastq.gz", 'w')

    with open(fq, 'r') as f:
        lines = []
        for line in f:
            lines.append(line.rstrip())
            if len(lines) == 4:
                record = process(lines)

                if record["sequence"][4:7] == barcodes[0]:
                    write_record(record, o1)

                if record["sequence"][4:7] == barcodes[1]:
                    write_record(record, o2)

                # sys.stderr.write("Record: %s\n" % (str(record)))
                lines = []
    o1.close()
    o2.close()

if __name__ == '__main__':
    print('{}\tStart!'.format(datetime.datetime.now()))
    fq = "/Users/shansabri/Desktop/APJ1_CLIPmm9_lane6_.fq"
    # fq = "/Users/shansabri/Desktop/test.fq"
    barcodes = ["AGT", "CCC"]
    demux(fq, barcodes)
    print('{}\tFinish!'.format(datetime.datetime.now()))