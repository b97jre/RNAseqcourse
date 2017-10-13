#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J STARref
#SBATCH -n 1

# load the sratoolkit module
module load sratoolkit/2.8.2

# run the sratoolkit fastq-dump.2.8.2 on the SRRfiles that you want to download
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199305
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199306
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199307
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199308
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199309
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199310
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199311
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199312
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199313
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199314
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199315
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199316
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199317
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199318
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199319
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199320
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199321
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199322
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199323
fastq-dump.2.8.2  -X 2000000 --split-3 SRR5199324


