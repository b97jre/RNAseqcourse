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

Reference guided assembly using Cufflinks
=========================================
Assembly of genes and isoforms using an assembly usig cufflinks is a two step process. 
First you map all the reads from your experiment to the reference sequence then you run another 
step where you use the mapped reads to idenitify potenital genes including introns and exons.  


Mapping short reads to a reference
----------------------------------

Here we will map the reads to the hg19 reference genome using a popular RNA-seq 
aligner, STAR. There are many many features that can be tweaked using STAR. We will use 
one flag **--outSAMstrandField intronMotif** which is needed for non-stranded RNAs to be handled 
properly by cufflinks.  

To load the STAR package on Uppmax, execute::

     module load star
     module load samtools

Now you can map the reads from one of the samples (or several; it's up to you 
which ones(s)) to map using a command such as the one below. ::
  
  mkdir sample1
    
  STAR  --genomeDir /proj/b2013006/downloads/courses/RNAseqWorkshop/reference/hg19_Gencode14.overhang75  --readFilesIn sample1_RAB11FIP5_1.fastq sample1_RAB11FIP5_2.fastq --runThreadN 2 --outSAMstrandField intronMotif --outFileNamePrefix sample1
	
flags used are 

* ``--genomeDir /path/to/STARgenome`` is the path to a pre-rendered reference library that STAR uses to map reads to the genome. 

*  ``--readFilesIn /path/to/read1 [/path/to/read2 ]`` is where you should add your fastq files that you will map to the reference.

*  ``--runThreadN N`` is the number of threads that will be used by the program.

*  ``--outSAMstrandField intronMotif`` adds SAM information so that unstranded reads function with cufflinks 

*  ``--outFileNamePrefix outDir`` tells you were the output data ends up. 


  
This should run fairly quickly and put a file called ``Aligned.out.sam`` in 
the directory that you specified with ``--outFileNamePrefix``. 

Next step is to convert it to a BAM file and rename it to a more representable name. 
A good naming praxis is to name the file to correspond what you mapped. As an example if you mapped sample 12
you should rename the mapped SAM file to a file with the name ``sample12_RAB11FIP5.bam`` 
The renaming and BAMfile conversion can be done in one step. Then to view them you also have to sort the hits and index them: ::

  
  samtools view -bSh -o sampleNN_RAB11FIP5.bam /path/to/Aligned.out.sam
  samtools sort sampleNN_RAB11FIP5.bam  sampleNN_RAB11FIP5.sorted
  samtools index sampleNN_RAB11FIP5.sorted.bam


The sorted, indexed bam file can then be viewed in IGV. 


Reference guided assembly using Cufflinks
-----------------------------------------

Cufflinks can do reference based assembly, which means 
that it tries to discover transcripts, disregarding gene annotation (actually there
is an option to use it as well but we will ignore that for now), just based on the 
mapped reads to the genome. This functionality works even on our small files.

Try to do a reference guided assembly. This is done simply by running Cufflinks 
without feeding it a GTF file with the -G flag::

     module load cufflinks/2.2.1
     cufflinks -o /path/to/outputDirectory sampleNN_RAB11FIP5.sorted.bam

Substitute the appropriate names for the BAM file and the output directory. When 
Cufflinks has finished (which should hardly take any time at all), the output 
directory will contain a file called ``transcripts.gtf``. Rename the file to a 
name that reflects what created the GTF file.  You can import that in 
the usual way into IGV as a track.

Was Cufflinks able to assemble your alignments into something that makes sense?
 
Other alternatives for reference-based assembly include 
`Scripture <http://www.broadinstitute.org/software/scripture>`_, 
`iReckon <http://compbio.cs.toronto.edu/ireckon/>`_ and 
`SLIDE <https://sites.google.com/site/jingyijli/>`_. These may require some 
annotation as input but they can discover (and quantify) new isoforms. 




Quantification with Cufflinks
=============================

Cufflinks is a well-known software package for estimating gene and isoform 
abundances in a BAM or SAM files based on an annotation file in GTF format 
(or for doing reference guided assembly which we will do below). However, we run 
into problems here because of the small file size. Cufflinks needs a certain amount 
of data to be able to do its estimations so it will not work very well on our small 
alignment files. Therefore we have run it on the full alignment files on each sample 
and provided the merged results at ``/proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/isoform_fpkm_table.txt``
(isoform FPKM) and ``/proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/RNAseqfpkm_table.txt``.
If you are interested, you can check in the isoform FPKM file how much each isoform 
is deemed to have been expressed by Cufflinks, but we will not go into that today. 
(These results will be revisited tomorrow, in the differential expression lab exercise.)

For reference, the commands we used were of the form::

     # Only for reference, does not need to be executed during the exercise
     cufflinks -p 8 -G /proj/g2013179/private/nobackup/RNAseqWorkshop/reference/Homo_sapiens.GRCh37.71.withchr.clean.gtf -o cufflinks_out_137_1 accepted_hits_137_1.bam

The ``-G`` option points to an annotation file in GTF format for which to calculate
FPKM values. The input here is a BAM file which is just a binary version of a SAM file.  

Other options for doing abundance estimation are `RSEM <http://deweylab.biostat.wisc.edu/rsem/>`_ 
or the flexible `RPKMforgenes.py script <http://sandberg.cmb.ki.se/media/data/rnaseq/instructions-rpkmforgenes.html>`_.





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



