## iCLIP

A few scripts to demultiplex and process iCLIP data

## Code Example

```bash
iCLIP_demultiplex.py --directory SxaQSEQsYA010L1 --barcodes AAG ACT ATC AGA GCC GTT --edits 0 --n-before-bc 4 --n-after-bc 4 --verbose

iCLIP_preprocess.py --fq-directory ~/scratch/fastq --adapter TGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT --umi-length 11 --min-length 20 --fastqc

iCLIP_pipeline.py --fq-directory ~/scratch/fastq --genome-dir ~/GSNAP/mm9 --genome mm9 --splice-directory ~/Genes/mm9.splicesites.iit --memory 4 --num-threads 6 --run-time 12
```
