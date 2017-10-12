#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J STARref
#SBATCH -n 16

# load the blast module
module load htseq/0.9.1


# run the STAR with 4 CPU threads (cores)
htseq-count -f bam -r pos -t exon -i gene -s no samFiles/CT_1_S1_L001_R1_0012.fastq.sorted.bam samFiles/CT_2_S2_L001_R1_0012.fastq.sorted.bam samFiles/CT_3_S3_L001_R1_0012.fastq.sorted.bam samFiles/CT_7_S7_L001_R1_0012.fastq.sorted.bam samFiles/CT_8_S8_L001_R1_0012.fastq.sorted.bam samFiles/CT_9_S9_L001_R1_0012.fastq.sorted.bam samFiles/WT_10_S10_L001_R1_0012.fastq.sorted.bam samFiles/WT_11_S11_L001_R1_0012.fastq.sorted.bam samFiles/WT_12_S12_L001_R1_0012.fastq.sorted.bam samFiles/WT_4_S4_L001_R1_0012.fastq.sorted.bam samFiles/WT_5_S5_L001_R1_0012.fastq.sorted.bam samFiles/WT_6_S6_L001_R1_0012.fastq.sorted.bam   GCF_000340665.1_C.cajan_V1.0_genomic.gene.gff 1>countTable_1.tab.txt 2>countTable_1.stderr &
htseq-count -f bam -r pos -t exon -i gene -s no samFiles/SRR4341972.sorted.bam samFiles/SRR4341973.sorted.bam samFiles/SRR4341974.sorted.bam samFiles/SRR4341975.sorted.bam samFiles/SRR4341976.sorted.bam samFiles/SRR4341977.sorted.bam samFiles/SRR4341978.sorted.bam samFiles/SRR4341979.sorted.bam  GCF_000340665.1_C.cajan_V1.0_genomic.gene.gff 1>countTable_2.tab.txt 2>countTable_2.stderr &
htseq-count -f bam -r pos -t exon -i gene -s no samFiles/SRR5199306.sorted.bam samFiles/SRR5199307.sorted.bam samFiles/SRR5199308.sorted.bam samFiles/SRR5199309.sorted.bam samFiles/SRR5199310.sorted.bam samFiles/SRR5199311.sorted.bam samFiles/SRR5199312.sorted.bam samFiles/SRR5199313.sorted.bam samFiles/SRR5199314.sorted.bam samFiles/SRR5199304.sorted.bam samFiles/SRR5199305.sorted.bam  GCF_000340665.1_C.cajan_V1.0_genomic.gene.gff 1>countTable_3.tab.txt 2>countTable_3.stderr &

wait
 
