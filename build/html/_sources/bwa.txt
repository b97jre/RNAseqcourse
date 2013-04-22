============================================
Mapping reads back to the assembly
============================================

Overview
======================

There are many different mappers available to map your reads back to the
assemblies. Usually they result in a SAM or BAM file
(http://genome.sph.umich.edu/wiki/SAM). Those are formats that contains the
alignment information, where BAM is the binary version of the plain text SAM
format. In this tutorial we will be using BWA
(http://bio-bwa.sourceforge.net/bwa.shtml).


The SAM/BAM file can afterwards be processed with Picard
(http://picard.sourceforge.net/) to remove duplicate reads. Those are likely to
be reads that come from a PCR duplicate (http://www.biostars.org/p/15818/).


BEDTools (http://code.google.com/p/bedtools/) can then be used to retrieve
coverage statistics.


Again, there is a script available that does it all at once. Read it and try to
understand what happens::
    
    less `which map-bwa-markduplicates`
    less /bubo/home/h16/inod/glob/github/metassemble/scripts/map/map-bwa-markduplicates.sh
    map-bwa-markduplicates -h

Mapping reads with bwa
======================
Take an assembly and try to map the reads back using bwa. Do this on an
interactive node again::

    # Create a new directory and link files
    mkdir -p ~/glob/asm-workshop/bwa
    cd ~/glob/asm-workshop/bwa
    ln -s ../velvet/out_21/contigs.fa contigs.fa
    ln -s ../sickle/pair1.fastq pair1.fastq
    ln -s ../sickle/pair1.fastq pair2.fastq

    # Run the everything in one go script. 
    map-bwa-markduplicates -t 8 -c pair1.fastq pair2.fastq pair contigs.fa contigs map


Some general statistics from the SAM/BAM file
=============================================
A good way to assess the quality of an assembly is by mapping the reads back
and determining what proportion was mapped. Use for instance::
    
    # Mapped reads only
    samtools view -c -F 4 map/contigs_pair-smds.bam
     
    # Unmapped reads only
    samtools view -c -f 4 map/contigs_pair-smds.bam

From: http://left.subtree.org/2012/04/13/counting-the-number-of-reads-in-a-bam-file/

Coverage information from BEDTools
=============================================
Look at the output from BEDTools::

    less map/contigs_pair-smds.coverage


Viewing the BAM file with Tablet
================================
Lastly one can look at the BAM file in tablet. If you have X-forwarding enabled
(``ssh -X inod@kalkyl1.uppmax.uu.se``) simply running::
    
    tablet

should do the trick. Please use an interactive node again.

(Optional) Run the mapping on all your contigs
===============================================
Run the mapping on all the assemblies and compare the number of mapped pairs in
each.
