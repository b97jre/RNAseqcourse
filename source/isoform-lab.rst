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

Mapping short reads to reference
=================================

Here we will map the reads to the hg19 reference genome using a popular RNA-seq 
aligner, STAR. There are many many features that can be tweaked using STAR. We will only use 

To load the TopHat package on Uppmax, execute::

     module load tophat/2.0.4
     module load samtools

Now you can map the reads from one of the samples (or several; it's up to you 
which ones(s)) to map using a command such as the one below. More important is that you understand (look at the 
the reason for using the flags you are using.  Note that you can use gzip-compressed 
files directly with Tophat (not that it is really needed for these tiny files; 
we just show it because it might be good to know.)::

     tophat --phred64-quals -p 8 -o my_test /sw/data/uppnex/reference/Homo_sapiens/hg19/program_files/bowtie2/concat /proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/sub_fastq/sample12_RAB11FIP5_1.fastq.gz /proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/sub_fastq/sample12_RAB11FIP5_2.fastq.gz

Since you should be on a full computing node with 4 processors, the flag ``-p 8`` 
will make the program run faster by using 8 threads. 

``--phred64-quals`` tells TopHat that these FASTQ files have their quality scores 
encoded in Illumina's older format called phred64, which is different from the standard Sanger (or phred33) format. 

``-o`` is the name of the output directory where the alignments will be written. 
You should put something that you can recognize as your own, like your name.

The long string ending in ``concat`` is the location of a Bowtie2 genome reference 
index, and the other options are file names of the two FASTQ files corresponding to 
the reads (each mate, read 1 and read 2, is in its own file.)

This should run fairly quickly and put a file called ``accepted_hits.bam`` in 
the directory that you specified with ``-o``. You should probably rename this 
file to something like ``sample12_RAB11FIP5.bam`` (if you mapped sample 12) to 
keep track of what the file represents. Now you can visualize the read "pile-up" 
in your genome browser. Before the visualization, you need to build a so-called 
BAM file index::

     samtools index sample12_RAB11FIP5.bam

(substitute the actual name of your file if it is different from ``sample12_RAB11FIP5.bam``.) 
Then, in IGV, just choose File > Load from File ... and select the BAM file. 
(If you are running IGV locally and did the mapping on Uppmax, you will have to 
download the BAM file and the corresponding index file (.bai) from your work folder 
to your own computer first.) When the track has loaded (actually several tracks) you 
will see the coverage and a representation of the actual reads themselves. How do 
they relate to the three unknown peptides?

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

Reference guided assembly using Cufflinks
=========================================

As mentioned above, Cufflinks can also do reference based assembly, which means 
that it tries to discover transcripts, disregarding gene annotation (actually there
is an option to use it as well but we will ignore that for now), just based on the 
mappings to the genome. This functionality works even on our small files.

Try to do a reference guided assembly. This is done simply by running Cufflinks 
without feeding it a GTF file with the -G flag::

     module load cufflinks/2.0.2
     cufflinks -o my_cuff_denovo_sample12 sample12_RAB11FIP5.bam

Substitute the appropriate names for the BAM file and the output directory. When 
Cufflinks has finished (which should hardly take any time at all), the output 
directory will contain a file called ``transcripts.gtf``. You can import that in 
the usual way into IGV (perhaps after renaming it into something less anonymous 
than ``transcripts.gtf``) as a track.

Was Cufflinks able to assemble your alignments into something that makes sense?

Other alternatives for reference-based assembly include 
`Scripture <http://www.broadinstitute.org/software/scripture>`_, 
`iReckon <http://compbio.cs.toronto.edu/ireckon/>`_ and 
`SLIDE <https://sites.google.com/site/jingyijli/>`_. These may require some 
annotation as input but they can discover (and quantify) new isoforms. 


