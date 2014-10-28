====================================
Identifying isoforms in RNA-seq data
====================================

One of the advantages of RNA-seq data compared to microarrays is that you get 
isoform-level information 'for free' in addition to gene-level information. 
In this exercise, we are going to look at some ways to identify and visualize isoforms.

As you have hopefully read on the introduction page, we will use RNA-seq and quantitative 
mass-spectrometry (MS) data sets from the A431 cell line. These data were measured by a 
group at SciLifeLab Stockholm in a 'proteogenomics' study where the aim was to discover 
new genes or gene variants by deep proteomic profiling using an MS method, and mapping 
the obtained peptides back to the genome. 
The RNA-seq data was obtained to see if there was RNA-level support for the predicted novel 
peptides and transcript variants. We will look at one loci that were flagged by the research 
group as being interesting, and see what the RNA-seq data look like for that gene.


Strategies for using the RNA-seq data
=====================================

There are different ways to find out how the RNA-seq data shows the RAB11FIP5 gene to 
be expressed. Roughly speaking, we can talk about three different strategies:

- Mapping the reads to a reference genome and quantifying the gene and isoform FPKMs

- Mapping the reads to a reference genome and assembling transcripts based on the mapping (reference guided assembly)

- Assembling the reads *de novo*, ignoring the reference genome for now

In order to make these steps computationally feasible during the lab, we have extracted 
only those sequences that mapped to the RAB11FIP5 gene in each sample. These "sub-FASTQ" 
files can be found in ``/proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/sub_fastq``.


To do the reference guided assembly yourself go to `Reference guided assembly using Cufflinks 
<https://export.uppmax.uu.se/b2013006/courses/RNAseq201410/build/html/courseSource/isoform-lab.html>`_. 
This link also contains information on how to quantify already annotated genes and isoforms.

To do the *de novo* assembly yourself go to `Isoform detection using RNA seq *de novo* Assembly 
<https://export.uppmax.uu.se/b2013006/courses/RNAseq201410/build/html/courseSource/isoform-denovo.html>`_.


============================
Visualising the isoform data
============================

For the gene **RAB11FIP5** you will now hopefully have generated your own data that you can look at. 
If everything worked you will now have 
* One BAM file with short reads mapped to the genome 
* One GTF file  containig differnt isoform of **RAB11FIP5** based on the short reads that mapped.
* One BAM file with **Trinity** assembled transcripts mapped to the genome

Since it takes time to generate all data we have already created other files that you can also download and view in your browser.

Importing reference based isoform info to the **RAB11FIP5** gene
================================================================




Importing de novo assembled transcripts mapped to the **RAB11FIP5** gene
========================================================================



Importing reference based isoform info to the genome
====================================================



Importing the peptide track to the **RAB11FIP5** gene or the genome                                                           
===================================================================
As mentioned above, we will compare some identified peptides from a mass-spectrometry 
experiment with RNA-seq data from the same cell line. Let's start by importing the track 
with identified peptides from the MS experiment. 



Importing the Pac bio reads mapped to the genome                                                         
================================================
Unfortunately there are no pacbio reads that mapped to the **RAB11FIP5** gene but if you 


If you are running a genome browser on Uppmax, you can find the file in the download area 

``/proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/RNAseqhuman_A431_global-TDA-FDR1pc_green-known_red-novel.bed`` 

If you are running locally, you can download the file from the 

`download area <https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/human_A431_global-TDA-FDR1pc_green-known_red-novel.bed>`_.

Then, in IGV, select File > Load from File ... and navigate to the BED file (on 
Uppmax according to above or to the location where you downloaded it locally). From 
the BED file extension, IGV will automatically know to color the track according to 
peptide status (green for annotated peptides, red for novel peptides.)




As mentioned above, we will compare some identified peptides from a mass-spectrometry 
experiment with RNA-seq data from the same cell line. Let's start by importing the track 
with identified peptides from the MS experiment. 

If you are running a genome browser on Uppmax, you can find the file in the download area 

``/proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/RNAseqhuman_A431_global-TDA-FDR1pc_green-known_red-novel.bed`` 

If you are running locally, you can download the file from the 

`download area <https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/human_A431_global-TDA-FDR1pc_green-known_red-novel.bed>`_.

Then, in IGV, select File > Load from File ... and navigate to the BED file (on 
Uppmax according to above or to the location where you downloaded it locally). From 
the BED file extension, IGV will automatically know to color the track according to 
peptide status (green for annotated peptides, red for novel peptides.)


Going to a locus
================

We will look at the RAB11FIP5 gene, which was highlighted in the proteomics experiments 
as having a couple of unannotated peptides. In IGV, click the textbox which has the word 
Go on its right hand side. Type RAB11FIP5 in it and press Enter.

Look at what you see in the display. What kind of peptides (previously known or novel) 
have been identified and how do they correspond to existing annotation?































Visualising the different 

Importing the peptide track                                                          
===========================

As mentioned above, we will compare some identified peptides from a mass-spectrometry 
experiment with RNA-seq data from the same cell line. Let's start by importing the track 
with identified peptides from the MS experiment. 

If you are running a genome browser on Uppmax, you can find the file in the download area 

``/proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/RNAseqhuman_A431_global-TDA-FDR1pc_green-known_red-novel.bed`` 

If you are running locally, you can download the file from the 

`download area <https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/human_A431_global-TDA-FDR1pc_green-known_red-novel.bed>`_.

Then, in IGV, select File > Load from File ... and navigate to the BED file (on 
Uppmax according to above or to the location where you downloaded it locally). From 
the BED file extension, IGV will automatically know to color the track according to 
peptide status (green for annotated peptides, red for novel peptides.)

Going to a locus
================

We will look at the RAB11FIP5 gene, which was highlighted in the proteomics experiments 
as having a couple of unannotated peptides. In IGV, click the textbox which has the word 
Go on its right hand side. Type RAB11FIP5 in it and press Enter.

Look at what you see in the display. What kind of peptides (previously known or novel) 
have been identified and how do they correspond to existing annotation?



