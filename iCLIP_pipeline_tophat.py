#!/usr/bin/python

__author__  = 'Shan Sabri'
__email__   = 'ShanASabri@gmail.com'
__date__    = '2016/10/03'

import datetime, argparse, os, sys

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
    parser = argparse.ArgumentParser(description='Process iCLIP fastq files.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--fq-directory', required=True, type=is_dir, default=os.getcwd(),
                        help='fastq directory (default: current working directory)')
    required.add_argument('-r', '--ref-directory', required=True, default=os.getcwd(),
                          help='Tophat2 reference/genome directory (default: current working directory)')
    required.add_argument('-gtf', required=True, help='GTF file path')

    parser.add_argument('-s', default='mm9', help='species flag for clipper (default: %(default)s)')
    parser.add_argument('--memory', type=int, default=4,
                        help='memory in gigabytes needed for alignment, per thread (default: %(default)s)')
    parser.add_argument('--num-threads', type=int, default=6,
                        help='the number of threads to use for alignment (default: %(default)s)')
    parser.add_argument('--run-time', type=int, default=24,
                        help='node compute time in hours needed for Hoffman2 (default: %(default)s hours)')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    return args

###########---------------------------------------###########

def init_dir(dirname):
    print '{}\tCreating necessary directories'.format(datetime.datetime.now()-start)
    base = os.path.dirname(dirname)
    if not os.path.exists(os.path.join(base, 'BAMs')): os.makedirs(os.path.join(base, 'BAMs'))
    if not os.path.exists(os.path.join(base, 'peaks')): os.makedirs(os.path.join(base, 'peaks'))

    return True

###########---------------------------------------###########

SCRIPT_TEMPLATE = """\
#!/bin/bash
source ~/.bash_profile
#$ -cwd
#$ -o {log}
#$ -e {log}
#$ -m n
#$ -l h_data={args.memory}G,h_rt={args.run_time}:00:00
#$ -pe shared {args.num_threads}

# ALIGN
{tophat2} \
    {tophat2_args} \
    -G {args.gtf} \
    -o {out} \
    {args.ref_directory} \
    {f} \
    > {aligned_log}.align.log

mv {out}/accepted_hits.bam {tmp}
"""

###########---------------------------------------###########

if __name__ == '__main__':
    args = parse_user_args()
    print '{}\tStart!'.format(datetime.datetime.now())

    init_dir(args.fq_directory)

    fastqs = [os.path.join(args.fq_directory, f) for f in os.listdir(args.fq_directory) if '.uniq.' in f]
    log = os.path.join(os.path.dirname(args.fq_directory), 'logs')
    bams = os.path.join(os.path.dirname(args.fq_directory), 'BAMs')
    samtools = os.popen('which samtools').readline().strip()
    tophat2 = os.popen('which tophat2').readline().strip()
    clipper = os.popen('which clipper').readline().strip()

    tophat2_args = '--num-threads {} --library-type fr-unstranded -N 0 -a 4 -x 1 -g 1 --no-coverage-search'.format(args.num_threads)
    clipper_args = '--premRNA --bonferroni'

    for f in fastqs:

        print '{}\tGenerating pipeline for {}'.format(datetime.datetime.now() - start, os.path.basename(f))
        clean_up = []
        tmp = f.replace(".fq.gz", ".bam")
        tmp = tmp.replace("fastq", "BAMs")
        aligned_log = tmp.replace('BAMs', 'logs')
        out = aligned_log.replace(".trimmed.uniq.bam", "")
        aligned_bam = f.replace(".fq.gz", ".bam")
        aligned_bam = aligned_bam.replace("fastq", "BAMs")
        aligned_sorted_bam = os.path.basename(f.replace(".fq.gz", ".sort"))
        flagstats = os.path.basename(f.replace(".fq.gz", ".stats.txt"))
        peaks = f.replace('.fq.gz', '.bed')
        peaks = peaks.replace('fastq', 'peaks')
        clean_up.extend((tmp, aligned_bam))
        clean_up = ' '.join(clean_up)

        script = SCRIPT_TEMPLATE.format(**globals())

        id = os.path.basename(f).split('.')[0]

        script_name = "iCLIP_%s" % id

        # print script

        with open(script_name, 'w') as f:
            f.write(script)
        # os.system("qsub %s" % script_name)
        # os.system("rm %s" % script_name)

    print '{}\tFinish!'.format(datetime.datetime.now())