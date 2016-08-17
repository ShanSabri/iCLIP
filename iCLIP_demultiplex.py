#!/usr/bin/python

__author__  = 'Shan Sabri'
__email__   = 'ShanASabri@gmail.com'
__date__    = '2016/08/13'
__version__ = '1.1'

###########---------------------------------------###########
#
# CHANGE LOG
# 2016-08-15    [1] added `--n-before-bc` and `--n-after-bc` to adjust barcode location
#               [2] kept UMI on read to properly collapse in preprocessing (iCLIP_preprocessing.py)
#
###########---------------------------------------###########

import datetime, argparse, shutil, glob, editdistance, gzip, os, sys

start = datetime.datetime.now()

###########---------------------------------------###########

def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = '{0} is not a directory'.format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

###########---------------------------------------###########

def parse_user_args():
    parser = argparse.ArgumentParser(description='Demultiplex iCLIP QSEQ files.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-d', '--directory', required=True, type=is_dir, default=os.getcwd(),
                        help='qseq directory (default: current working directory)')
    required.add_argument('-b', '--barcodes', nargs='+', type=str, required=True,
                        help='a space-seperated list of barcodes to preform demultiplexing')
    parser.add_argument('-e', '--edits', type=int, default=0,
                        help='number of sequence edits to allow within the barcode (default: %(default)s)')
    parser.add_argument('--n-before-bc', type=int, default=4,
                        help='number of randomers before the barcode\n(default: %(default)s, nnnnXXX)')
    parser.add_argument('--n-after-bc', type=int, default=4,
                        help='number of randomers after the barcode\n(default: %(default)s, XXXnnnn)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='increase output verbosity')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    return args

###########---------------------------------------###########

def init_dir(dirname):
    print '{}\tCreating necessary directories'.format(datetime.datetime.now()-start)
    base = os.path.dirname(dirname)
    if not os.path.exists(os.path.join(base, 'fastq')): os.makedirs(os.path.join(base, 'fastq'))
    if not os.path.exists(os.path.join(base, 'tmp')): os.makedirs(os.path.join(base, 'tmp'))

    return True

###########---------------------------------------###########

def cat_qseqs(dirname):
    print '{}\tConcatenating QSEQ files'.format(datetime.datetime.now() - start)
    prefix = set(['_'.join(os.path.splitext(f)[0].split('_')[0:3]) for f in os.listdir(dirname) if f.endswith('.txt.gz')])
    for p in prefix:
        files = [os.path.join(dirname, f) for f in os.listdir(dirname) if f.startswith(p)]
        with open(os.path.join(os.path.dirname(dirname), 'tmp', p) + '.qseq.txt.gz', 'wb') as wfp:
            for fn in files:
                with open(fn, 'rb') as rfp:
                    shutil.copyfileobj(rfp, wfp)

    return True

###########---------------------------------------###########

def demux(args):
    cat_qseqs(args.directory)

    print '{}\tDemultiplexing QSEQ files'.format(datetime.datetime.now() - start)
    out = {}
    for t in args.barcodes: out[t] = gzip.open(os.path.join(os.path.dirname(args.directory), 'fastq', t) + '.fq.gz', 'w')

    largest = sorted((os.path.getsize(s), s) for s in glob.glob(os.path.join(os.path.dirname(args.directory), 'tmp', '*')))[-1][1]
    with gzip.open(largest) as f:
        for line, f in enumerate(f, start=1):
            f = f.strip().split('\t')
            obs_idx = f[8][args.n_before_bc:(args.n_before_bc + 3)]
            obs_read = f[8] # [(args.n_before_bc + 3 + args.n_after_bc):]
            obs_read_qual = f[9] # [(args.n_before_bc + 3 + args.n_after_bc):]
            obs_read = obs_read.rstrip('.')
            obs_read_qual = obs_read_qual[:len(obs_read)]
            if len(obs_read) < 1 : continue
            obs_read = obs_read.replace(".", "N")
            obs_read_id = '@' + ':'.join(f[0:8]) + ' length:' + str(len(obs_read))

            def match_keys(_true_idx):
                return editdistance.eval(_true_idx, obs_idx)

            true_idx = min(out.keys(), key=match_keys)

            if editdistance.eval(true_idx, obs_idx) > args.edits:
                line=line+1
                continue

            if args.verbose:
                print '{}\t\tProcessing line {}: {} - {} to {}'.format(datetime.datetime.now() - start, line, ':'.join(f[0:8]), obs_idx, true_idx)

            out_file = out[true_idx]
            out_file.write('{0}\n{1}\n+\n{2}\n'.format(obs_read_id, obs_read, obs_read_qual))

    return True

###########---------------------------------------###########

def clean_up():
    dir = os.path.join(os.path.dirname(args.directory), 'tmp')
    shutil.rmtree(dir)

###########---------------------------------------###########

if __name__ == '__main__':
    args = parse_user_args()
    print '{}\tStart!'.format(datetime.datetime.now())
    init_dir(args.directory)
    demux(args)
    clean_up()
    print '{}\tFinish!'.format(datetime.datetime.now())