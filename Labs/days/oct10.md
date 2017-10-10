---
title: "10 oct workshop line of code"
author: "Johan"
date: "10/10/2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Cmds run on terminal



Make directory:

    mkdir RNASeqCourse2017

go to folder:  

    cd RNASeqCourse2017/



# Run fastqc 

Copy files from RNAseq course folder:

    cp /home/jreimegard/RNAseqCourse/QC/S15_S15_L001_R* .

  
Load the fastqc module:
    
    module load fastqc/0.11.5
  
run fastqc:
  
    fastqc  fastqc_output S15_S15_L001_R1_001.fastq.gz S15_S15_L001_R2_001.fastq.gz  

Wiew fastq file

In linux or mac os x or windows:  
    
    mkdir ~/public_html/
    #Copy files to folder ( linux and mac os x )
    cp *.html ~/public_html/

View the files in the browser:
  http://hpc.ilri.cgiar.org/~jreimegard/S15_S15_L001_R2_001_fastqc.html
  
  
  
  
# Removing adapter sequences  
If you want remove the adapter sequences using cutadapt:

    module load cutadapt/v1.8.1 
  
    cutadapt -a TACGGAGGATCCAAGCGTTATCCGGAATCATTGGGTTTAAAGGGTCCGTA -o S15_S15_L001_R1_001.cutadapt.fastq.gz S15_S15_L001_R1_001.fastq.gz 


# Map reads 

create folder:

    mkdir mapping
    cd mapping
    
    
copy fastqfiles to folder:

    cp  /home/jreimegard/RNAseqCourse/isoform/referenceBased/data/sample1_RAB11FIP5_* .
    
Link (soft) folder with hisat genome reference

    ln -s /home/jreimegard/RNAseqCourse/hg19_hisat2

Load hisat module and run hisat2

    mkdir samfiles
    module load hisat2/2.0.5
    hisat2 -p 1 --dta-cufflinks -x hg19_hisat2/genome -1 sample1_RAB11FIP5_1.fastq -2 sample1_RAB11FIP5_2.fastq -S samfiles/sample1.RAB11FIP5.sam &> samfiles/sample1.RAB11FIP5.sam.info
    
    
Load samtools and convert from sam to sorted bam and index it:

    cd samfiles
    module load samtools
    samtools view -bS -o sample1.RAB11FIP5.bam sample1.RAB11FIP5.sam 
    samtools sort sample1.RAB11FIP5.bam > sample1.RAB11FIP5.sorted.bam
    samtools index sample1.RAB11FIP5.sorted.bam 
    
    
Copy sorted bam file so it can be accessed via html 

    cp sample1.RAB11FIP5.sorted.ba* ~/public_html/
    
    
    


