#!/usr/bin/python

__author__  = 'Shan Sabri'
__email__   = 'ShanASabri@gmail.com'
__date__    = '2016/08/15'

import datetime, argparse, os, sys, gzip
from collections import Counter

start = datetime.datetime.now()

###########---------------------------------------###########
#
# CHANGE LOG
# 2016-08-22    [1] uniq read WITH umi, then write out without UMI
#               [2] write out fastq with arbitrary quality score
#               [3] added UMI to fastq header
#
###########---------------------------------------###########

###########---------------------------------------###########

def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = '{0} is not a directory'.format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

###########---------------------------------------###########

def parse_user_args():
    parser = argparse.ArgumentParser(description='Process iCLIP fastq files.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--fq-directory', required=True, type=is_dir, default=os.getcwd(),
                        help='fastq directory (default: current working directory)')
    required.add_argument('-a', '--adapter', required=True, type=str,
                          help='adapter to be clipped')
    parser.add_argument('--umi-length', type=int, default=11,
                        help='length of the UMI sequence (default: %(default)s)')
    parser.add_argument('-m', '--min-length', type=int, default=20,
                        help='after adapter trimming, toss sequences less than a specified length (default: %(default)s)')
    parser.add_argument('--fastqc', action='store_true',
                        help='run FastQC before and after trimming')
    parser.add_argument('--clean-up', action='store_true',
                        help='delete intermediate files (raw fastqs and trimmed fastqs, recommended if storage space is an issue)')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    return args

###########---------------------------------------###########

def init_dir(args):
    print('{}\tCreating necessary directories'.format(datetime.datetime.now()-start))
    base = os.path.dirname(args.fq_directory)
    if not os.path.exists(os.path.join(base, 'logs')): os.makedirs(os.path.join(base, 'logs'))
    if args.fastqc:
        if not os.path.exists(os.path.join(base, 'fastqc')): os.makedirs(os.path.join(base, 'fastqc'))

    return True

###########---------------------------------------###########

def run_fastqc(f):
    print('{}\tRunning FastQC on {}'.format(datetime.datetime.now() - start, os.path.basename(f)))
    fastqc = os.popen('which fastqc').readline().strip()
    if not fastqc:
        fastqc = '/Applications/FastQC.app/Contents/MacOS/fastqc'
    outdir = os.path.join(os.path.dirname(os.path.dirname(f)), 'fastqc')
    cmd = [fastqc, '--outdir ' + outdir, '--quiet', f]
    # print " ".join(cmd)
    output = os.popen(" ".join(cmd)).read()

    return True

###########---------------------------------------###########

def cutadapt(f, args):
    print('{}\tClipping adapter sequence from {}'.format(datetime.datetime.now() - start, os.path.basename(f)))
    cutadapt = os.popen('which cutadapt').readline().strip()
    output = f[:-6] + '.trimmed.fq.gz'
    log = os.path.join(os.path.dirname(args.fq_directory), 'logs', os.path.basename(f).split('.')[0] + '.trim.log')
    cmd = [cutadapt, "-a " + args.adapter,
           "--length-tag 'length='",
           "--minimum-length " + str(args.min_length),
           "-o " + output, f]
    # print " ".join(cmd)
    output_log = os.popen(" ".join(cmd)).read()
    with open(log, 'w') as out: out.write(output_log)

    return output

###########---------------------------------------###########

def next_block(fastq_file):
    l1 = fastq_file.readline()
    l2 = fastq_file.readline()
    l3 = fastq_file.readline()
    l4 = fastq_file.readline()
    return (l1,l2,l3,l4)

def valid_block(b):
    return b[0] != ''

def uniq_fq(f, args):
    print('{}\tRemoving duplicate reads from {}'.format(datetime.datetime.now() - start,  os.path.basename(f)))

    fastq_file = gzip.open(f)
    sequences = Counter()
    block = next_block(fastq_file)
    while valid_block(block):
        seq = block[1]
        if seq in sequences:
            sequences[seq] += 1
        else:
            sequences[seq] = 1
        block = next_block(fastq_file)
    fastq_file.close()

    output = f[:-6] + '.uniq.fq.gz'
    too_short_after_umi_cut = 0
    with gzip.open(output, 'wb') as out:
        for n, (k, v) in enumerate(sequences.most_common(), start=1):
            if len(k[args.umi_length:]) - 1 >= args.min_length:
                qual = 'D'*(len(k[args.umi_length:])-1)
                out.write('@Sequence_{}_{}_with_{}_occurrences\n{}\n+\n{}\n'.format(str(n), k[:args.umi_length], str(v),
                                                                                 k[args.umi_length:].strip(), qual))
            else:
                too_short_after_umi_cut += 1
                continue

    print('{}\tFound {} unique sequences in {} (total={})'.format(datetime.datetime.now() - start,
                                                                      len(sequences.keys()),  os.path.basename(f), sum(sequences.values())))
    print('{}\tFound the most common unique sequence to be {}'.format(datetime.datetime.now() - start,
                                                                          " ".join('{} occurring {} times'.format(k.strip(), str(v))
                                                                                   for k, v in sequences.most_common(1))))
    print('{}\t{} sequences failed to write because the length without the UMI is less than {}'.format(datetime.datetime.now() - start,
                                                                                                           too_short_after_umi_cut,
                                                                                                           args.min_length))
    return output

###########---------------------------------------###########

def process(args):
    fastqs = [os.path.join(args.fq_directory, f) for f in os.listdir(args.fq_directory) if f.endswith('.fq.gz')]
    os.system('source ~/.bash_profile')
    for f in fastqs:
        if args.fastqc: run_fastqc(f)
        clipped_fq = cutadapt(f, args)
        if args.fastqc: run_fastqc(clipped_fq)
        uniqued_fq = uniq_fq(clipped_fq, args)
        if args.fastqc: run_fastqc(uniqued_fq)

    return True

###########---------------------------------------###########

def clean_up(args):
    print('{}\tCleaning up!'.format(datetime.datetime.now() - start))
    intermediate_fq = [os.path.join(args.fq_directory, f) for f in os.listdir(args.fq_directory) if not '.uniq.' in f]
    fastqc_dir = os.path.join(os.path.dirname(args.fq_directory), 'fastqc')
    fastqc_zip = [os.path.join(fastqc_dir, f) for f in os.listdir(fastqc_dir) if f.endswith('.zip')]
    delete = intermediate_fq + fastqc_zip

    for f in delete: os.remove(f)

###########---------------------------------------###########

if __name__ == '__main__':
    args = parse_user_args()
    print('{}\tStart!'.format(datetime.datetime.now()))
    init_dir(args)
    process(args)
    if args.clean_up: clean_up(args)
    print('{}\tFinish!'.format(datetime.datetime.now()))
