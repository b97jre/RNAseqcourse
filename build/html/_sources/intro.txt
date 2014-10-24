===================
General information
===================
For this lab, we are going to work on UPPMAX, in the project ``g2014046``, which has been set up for the purposes of this course. You reach the home directory of the project by changing working directory to ``/proj/g2014046/nobackup/private/workshop/``::

    cd /proj/g2014046/nobackup/private/workshop/

The files we are going to use in the exercises are in ``/proj/g2014046/webexport/files/RNAseqWorkshop/download/``. This document (the web-based workshop documentation) is hosted under ``/proj/g2014046/webexport/RNAseq-workshop/build``. 

Some of the files that we will use have been pre-constructed by us ahead of the lab because it would take too much time to do it during the lab exercise. In those cases, we have indicated the commands that we used to generate the files, so that you should be able to reproduce the results on your own.


Short summary of the exercises
==============================

In the isoform exercise, and part of the differential expression exercise, we are going to look at some RNA-seq (and mass spectrometry) data from the A431 cell line. It is an epidermoid carcinoma cell line which is often used to study cancer and the cell cycle, and as a sort of positive control of epidermal growth factor receptor (EGFR) expression. A431 cells express very high levels of EGFR, in contrast to normal human fibroblasts. In the experiment we are looking at today, A431 cells were treated with gefinitib, which is an EGFR inhibitor, and is used (under the trade name Iressa) as a drug to treat cancers with mutated and overactive EGFR. In the experiment, RNA was extracted at four time points: before the gefinitib treatment (t=0), and two, six and twenty-four hours after treatment (t=2, t=6, t=24, respectively), and sequenced using an Illumina HiSeq instrument in triplicates (thus there are 3x4=12 samples). In a parallel experiment, protein contents at those same time points were measured using mass spectrometry in duplicates (resulting in 8 samples). 

There are many relevant questions that could be asked based on these measurements. During the isoform exercise, we are going to look at some specific regions where the mass-spectrometry data indicated that novel exons or splice variants could be present at the protein level. We will use (part of) the RNA-seq data to examine if there is corresponding evidence on the mRNA level, and how different software tools could be used to detect and visualize novel gene variants. 
