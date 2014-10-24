=============================
Short summary of the datasets
=============================
For the tutorials we have set up we have two different datasets.
The first set are polyA enriched RNAsequenses from a A431 cell line.     
The second set are sRNA sequences from drosophila

Short summary of the polyA enriched dataset
===========================================

In the isoform exercise, and part of the differential expression exercise, 
we are going to look at some RNA-seq data from the A431 cell line. 
It is an epidermoid carcinoma cell line which is often used to study cancer
 and the cell cycle, and as a sort of positive control of epidermal growth factor
  receptor (EGFR) expression. A431 cells express very high levels of EGFR, in contrast
 to normal human fibroblasts. 
 
 In the experiment we are looking at today, A431 cells were treated with gefinitib, which is an EGFR inhibitor, 
 and is used (under the trade name Iressa) as a drug to treat cancers with mutated and overactive EGFR. 
 In the experiment, RNA was extracted at four time points: before the gefinitib treatment (t=0), and two, six
 and twenty-four hours after treatment (t=2, t=6, t=24, respectively), and sequenced using an Illumina
 HiSeq instrument in triplicates (thus there are 3x4=12 samples). 

We also have peptide data mapped to the genome sequence that can support the presence of exons in the sequence. 
 
This dataset or parts of it will be used in the QC lab, the isoform labs and the differential expression lab.
There are many relevant questions that could be asked based on these measurements. 
In the QC exercise we are going to examine if the RNA libraries that we work with are what we think they are or if 
there are some missannotations on the datasets.
In the isoform exercises we are going to look at some specific regions where the mass-spectrometry data 
indicated that novel exons or splice variants could be present at the protein level. We will use (part of) 
the RNA-seq data to examine if there is corresponding evidence on the mRNA level, 
and how different software tools could be used to detect novel gene variants. 

 
Short summary of the sRNA dataset
=================================
In this exercise we will analyze a few small RNA libraries, from Drosophila melanogaster (fruit fly) embryos
and two cell lines (KC167 cells derived from whole embryos, and ML-DmD32 cells derived from adult wing discs).
This is a subset of a much larger data set used to study microRNAs and other small RNAs in Drosophila
These data sets are described more in this paper: http://genome.cshlp.org/content/24/7/1236.full

Where to find the data 
======================
The data for each lab is can be downloaded using two different methods. 
