---
title: "STAR reference"
author: "Johan"
date: "10/10/2017"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Create a sbatch script to create a STAR reference


Download reference to your folder





Make directory for the star reference:

    mkdir referenceGenomeSTAR

Create a text file **runSTARreference.sbatch** with the following information:

    #!/usr/bin/env bash
    #SBATCH -p batch
    #SBATCH -J STARref
    #SBATCH -n 8

    # load the blast module
    module load star/2.5.3a

    # run the STAR with 4 CPU threads (cores)
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir pigeonPea_STAR --genomeFastaFiles GCF_000340665.1_C.cajan_V1.0_genomic.fna



Start the sbatch script job with:

    sbatch runSTARreference.sbatch