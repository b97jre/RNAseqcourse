================================================
Isoform detection using RNA seq de novo Assembly 
================================================

We are going to use one of the open source RNA *de novo* assemblers 
during this practical. It is called **Trinity**. Independent assessment 
of *de novo* assembly programs showed that Trinity was one of the best assemblers to use. 
It is also one of the programs that is being updated and does also have downstream analysis tools. 

.. image:: ../images/Compare.jpg

Figure taken from `Optimizing de novo transcriptome assembly from short-read RNA-Seq data: a comparative study 
<http://www.biomedcentral.com/1471-2105/12/S14/S2>`_.

A de novo  take your reads and turn them into *contigs*. For more details
on how **Trinity** work read the corresponding paper (`Trinity 
<http://www.nature.com/nbt/journal/v29/n7/full/nbt.1883.html>`_
). 


Running an assembly using **Trinity** for ~20 million of reads takes at least a day and more than 50 GB of RAM. In order 
to reduce the time for the **Trinity** to run during this course we will focus on reads that can be mapped to a small region on the human chromosome.  


Preparation
===========

If you are working on UPPMAX i suggest that you will do all your exercises from your glob folder. 


Make a new subdirectory and go there for this exercise.  ::


   #Go to your own glob folder mkdir deNovoAssembly  
   cd ~/glob/RNAseqWorkshop
   
   #Create a new folder where you will do this exercise
   mkdir deNovoAssembly  
   # go into that directory 
   cd deNovoAssembly
   
   
   
Files used during the exercise 
==============================
 
   
If you are on uppmax
--------------------

Copy all the files that you will need for this exercise from here. ::

    # This assumes that you now are in the deNovoAssembly folder 
    cp -r /proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/deNovo/data . 


If you are somewhere else
-------------------------
You can download all data using a webinterface from `here
<https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/deNovo/>`_ 


   
Programs used during the exercise 
=================================

When doing this course on UPPMAX all programs will be available to load as pre-installed modules. 
In order to be able to use these programs you need to load the modules before using them. 

Load all programs that you will need for trinity to work on uppmax this exercise. ::
 
    # load modules to make trinity work 
    module load bioinfo-tools 
    module load bowtie
    module load samtools
    module load trinity/2014-07-17 
    
    
Load all programs that you will need for STAR to work on UPPMAX this exercise. ::

    # load modules to make RNAseq aligner STAR work 
    
    module load star

   
If you are doing this exercise on somewhere else follow each program information on how to install it.
   

Assemble the reads into contigs 
===============================

Since **Trinity** is often being updated you should make sure you are using the latest version.
That means that the requirements and the command line to to run **Trinity** changes occasionally. 
You can find the basic usage info for **Trinity** `here
http://trinityrnaseq.sourceforge.net/#running_trinity
 and the latest typical command line to type `here
<http://trinityrnaseq.sourceforge.net/#running_trinity>`_ 

Go to the webpage and try out the code that is 

Good things to know about the data used in this lab. ::

    # You will use paired end data. 
    # The file that trinity refers to as left is the fastq files that ends with _1.fastq
    # The file that trinity refers to as right is the fastq files that ends with _2.fastq
 
    # The RNA seq data that we use in this exercise is not strand specific.
     
	 


To fully use the potential 
of the programs it is worthwhile to read the manual and use the correct flags. As 
an example both programs handle strand specific RNA that reduces the complexity of 
the algorithm and also produces better results.
or `Trinity
<http://trinityrnaseq.sourceforge.net/#running_trinity>`_
manuals.





Assessing the new assemblies
===========================

Now that the reads have been assembled into contigs you can map them back onto 
the human genome sequence to see how they were assembled. Note that in 
non-model organism this is not possible. If you would like to asses the assembly
of transcripts without a reference genome Trinity has a downstream analysis pipeline 
that is worth following `Trinotate
<http://trinityrnaseq.sourceforge.net/annotation/Trinotate.html>`_ . This is not something we will 
do in this course but if you have time over feel free to try it out. 
    
Start with mapping the trinity assembled transcripts to the human genome using STAR. 
Convert them to bam format, sort and index them using samtools::
  
  mkdir STARtrinityMapping
    
  STAR  --genomeDir /proj/g2014046/private/RNAseqWorkshop/reference/hg19_Gencode14.overhang75  --readFilesIn Trinity/Trinity.fasta --runThreadN 1 --outSAMstrandField intronMotif --outFileNamePrefix STARtrinityMapping/
  samtools view -bSh -o trinityTranscripts.bam STARtrinityMapping/Aligned.out.sam
  samtools sort trinityTranscripts.bam  trinityTranscripts.sorted
  samtools index trinityTranscripts.sorted.bam
	
Do the same procedure for the oases assembled transcripts::
	
  mkdir STARoasesMapping
  STAR  --genomeDir /proj/g2014046/private/RNAseqWorkshop/reference/hg19_Gencode14.overhang75  --readFilesIn oasesPipelineMerged/transcripts.fa --runThreadN 1 --outSAMstrandField intronMotif --outFileNamePrefix STARoasesMapping/
  samtools view -bSh -o oasesTranscripts.bam STARoasesMapping/Aligned.out.sam
  samtools sort oasesTranscripts.bam  oasesTranscripts.sorted
  samtools index oasesTranscripts.sorted.bam
	
    
When ready there should be two bam files that are sorted and indexed. These can now be viewed in the IGV 
or Savant genome browsers. In total there were 12 samples and you have now assembled one of those samples. 
If time permits do one more sample. If time is running out you can download and view all the 24 different samples. 
We have also merged the reads from all the 12 samples and used all the reads to create assembled transcripts.
All these files can be found `here 
<https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/deNovoFinishedFiles/AllBamFiles/>`_

Download a few of them and compare the differents states to see if you can identify different isoforms. How does the 
de novo assembled transcripts compare to the reference based isoform detection programs. 
    


**OPTIONAL**
I recomend to download the bamFiles and view them in a genome browser on your laptop.
The interactive genome view experience on UPPMAX, especially when loading many tracks, can 
be slow.This is done in two steps. ::

    #create a folder for all the bamfiles
    mkdir AllBamFiles 
    
    # move all the bamfiles into that folder 
    mv *.sorted.bam AllBamFiles
    
    #create a tar file with all the bamFiles so that you can download them to your laptop

    tar -cf AllBamFiles.tar AllBamFiles 
    
    #Use any sftp program of your choice to download the files from uppmax
    
    # If you are using shell you can open up a new terminal window and go to 
    # the place where you want to store your bamFiles
    
    cd $YOURLOCALPATH
    scp yourUppmaxName@milou.uppmax.uu.se:$WORKDIR/deNovo/AllBamFiles.tar . 
    
    
    
    
Now that you have all the bam files in with individual names try to view them in 
a genome brower, both IGV and Savant works fine. Here we will describe how to view them
in IGV but SAVANT has a nice feature of viewing paired end reads as arcs that IGV
misses. If you have time i recomend trying both of them out. 

First have a look on the  two bamfiles that contains the assemblies of all
reads from all twelve timepoints with the two different assemblers. They have the 
names ``RAB11FIP5_trinity.Trinity._hg_19_STAR.bam`` and 
``RAB11FIP5_oases.Oases._hg_19_STAR.bam``. ::

    #If you view your files on your laptop start IGV like this

    java -Xmx1500M -jar igv.jar
    
    # If you view your files on UPPMAX do according to UPPMAX
    
    
    #Load tracks in the IGV browser
    
    File->Load From File...
    	choose **oasesTranscripts.sorted.bam**
    	and **trinityTranscripts.sorted.bam**
    	
    # Load peptide sequences 	

    File->Load From File...
    	choose **human_A431_global-TDA-FDR1pc_green-known_red-novel.bed**
    	
    # Load your mapped reads from before  	

    File->Load From File...
    	choose **sample12_RAB11FIP5.bam**
    	
    # Load your own GTF file
    
    File->Load From File...
    	choose **transcripts.gtf** or what you have named it.
    	

**OPTIONAL**
There is also a possibility to view tracks that is publicly available. This is easy to 
do in IGV and adds some information in the region that we are looking into. ::
	
    	
    # Load different gene annotations files

    File->Load From Server...
    	choose Available Datasets ->Annotations -> Genes ->UCSC Genes
    

    # Load multiple alignments to other vertebrates

    File->Load From Server...
    	choose Available Datasets ->Annotations -> Comparative Genomics ->Phastcons (Vertebrate 46 way)
	

   # Load any of the other annotations that you think is interesting

    File->Load From Server...
    	choose Available Datasets ->..  -> ..  ->Up to you 
	



    	
Now have a look at the de novo assembled transcripts. Do they seem reasonable? Which 
regions on the de novo assembled transcripts do not correspond to your own .gtf 
file?  Which is the correct one? 

Now take a closer look at the region chr2:73,308,166-73,308,278. This corresponds 
to the regions where the RefSeq genes is annotated as intron but the *de novo* assembly
, the cufflinks gtf file and the peptide file suggest that the region is being transcribed 
and translated into peptides. When examining the *de novo* assembled contigs it seems
that none of the transcripts goes through the region. Is this real or could there 
be a shortcoming of the assembler or the sequencing platform? Unfortunately we do 
not have the answers to these questions but all the different methods add in to give 
more understanding in the complexity of isoform analysis and genome annotation.  
    	
    
    
    
    
    
	
	
	
   
     
	



	





