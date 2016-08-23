#!/usr/bin/python

__author__  = 'Shan Sabri'
__email__   = 'ShanASabri@gmail.com'
__date__    = '2016/08/16'

import datetime, argparse, os, sys

start = datetime.datetime.now()

###########---------------------------------------###########
#
# Hoffman2 usage example
#
# python iCLIP_pipeline.py
#     --fq-directory /u/home/s/ssabri/scratch/Martina_Roos/2016-08-12/fastq
#     --genome-dir /u/home/s/ssabri/project-ernst/ref_data/Mus_musculus/UCSC/mm9/GSNAP/mm9
#     --genome mm9 --splice-directory /u/home/s/ssabri/project-ernst/ref_data/Mus_musculus/UCSC/mm9/Annotation/Genes/mm9.splicesites.iit
#     --memory 4 --num-threads 6 --run-time 12
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
    required.add_argument('--genome-directory', required=True,
                          help='built genome directory')
    required.add_argument('--splice-directory', required=True,
                        help='splicing <STRING>.iit directory')
    parser.add_argument('-g', '--genome', type=str, default='mm9',
                        help='built genome database (default: %(default)s)')
    parser.add_argument('--memory', type=int, default=4,
                        help='memory in gigabytes needed for alignment, per thread (default: %(default)s)')
    parser.add_argument('--num-threads', type=int, default=8,
                        help='the number of threads to use for alignment (default: %(default)s)')
    parser.add_argument('--run-time', type=int, default=24,
                        help='node compute time in hours needed for Hoffman2 (default: %(default)s hours)')
    parser.add_argument('--gsnap-args', nargs='+', type=str,
                          help='a space-seperated list of arguments for GSNAP (default: -t 4 -N 1 --max-mismatches=0 -A sam --gunzip)')
    parser.add_argument('--clipper-args', nargs='+', type=str,
                        help='a space-seperated list of arguments for CLIPPER (default: --premRNA --bonferroni)')

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
{gsnap} \
    {gsnap_args} \
    -s {args.splice_directory} \
    -D {args.genome_directory} \
    -d {args.genome} \
    {f} \
    > {aligned_sam} \
    2> {aligned_log}.align.log


# CONVERT SAM TO BAM
{samtools} \
    view \
    -bS \
    {aligned_sam} \
    > {aligned_bam}

# SORT BAM
{samtools} \
    sort \
    {aligned_bam} \
    {bams}/{aligned_sorted_bam}

# FLAGSTATS
{samtools} \
    flagstat \
    {bams}/{aligned_sorted_bam}.bam \
    > {log}/{flagstats}

# INDEX
{samtools} \
    index \
    {bams}/{aligned_sorted_bam}.bam

# CLIPPER
{clipper} \
    -b {bams}/{aligned_sorted_bam}.bam \
    -s {args.genome} \
    {clipper_args} \
    --outfile {peaks}

# CLEAN UP
rm {clean_up}
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
    gsnap = os.popen('which gsnap').readline().strip()
    clipper = os.popen('which clipper').readline().strip()

    if args.gsnap_args:
        gsnap_args = ' '.join(args.gsnap_args)
    else:
        gsnap_args = '-t {} -N 1 -A sam --gunzip -B 5'.format(args.num_threads)

    if args.clipper_args:
        clipper_args = ' '.join(args.clipper_args)
    else:
        clipper_args = '--premRNA --bonferroni'

    for f in fastqs:
        print '{}\tGenerating pipeline for {}'.format(datetime.datetime.now() - start, os.path.basename(f))
        clean_up = []
        aligned_sam = f.replace(".fq.gz", ".sam")
        aligned_sam = aligned_sam.replace("fastq", "BAMs")
        aligned_log = aligned_sam.replace('BAMs', 'logs')
        aligned_bam = f.replace(".fq.gz", ".bam")
        aligned_bam = aligned_bam.replace("fastq", "BAMs")
        aligned_sorted_bam = os.path.basename(f.replace(".fq.gz", ".sort"))
        flagstats = os.path.basename(f.replace(".fq.gz", ".stats.txt"))
        peaks = f.replace('.fq.gz', '.bed')
        peaks = peaks.replace('fastq', 'peaks')
        clean_up.extend((aligned_sam, aligned_bam))
        clean_up = ' '.join(clean_up)

        script = SCRIPT_TEMPLATE.format(**globals())

        id = os.path.basename(f).split('.')[0]

        script_name = "iCLIP_%s" % id

        # print script

        with open(script_name, 'w') as f:
            f.write(script)
        os.system("qsub %s" % script_name)
        os.system("rm %s" % script_name)

    print '{}\tFinish!'.format(datetime.datetime.now())