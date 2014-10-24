=======================
Links and Files and FAQ
=======================


Working directory
=================

Directory where everyone should keep their data and work

``/proj/g2014046/nobackup/private/workshop``

Everyone should create their own folder with the same name as their uppmax username were they should do all the work during 
the workshop.

To create your own subdirectory under workshop do the following::

    cd /proj/g2014046/nobackup/private/workshop
    
    #create a folder with your uppmax user account namne
    # in my case it is johanr you should use your name
    
    MYNAME=johanr
    
    mkdir $MYNAME
    




Data directories
================

Directory were all data that will be used during the workshop is located:

RNA seq part (Exercises on tuesday and wednesday)

`Link to RNA-seq exercise download area 
<https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/>`_

Same files can be found at uppmax at the following directory

``/proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq``

Proteomics part (Exercises on thursday)

`Link to Proteomics exercise download area
<https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/Proteomics/>`_

Same files can be found at uppmax at the following directory

``/proj/g2014046/webexport/files/RNAseqWorkshop/download/Proteomics``


FAQ
===

**Link is not working on webpage**

Tell Johan that something is wrong and wait a few minutes and then update the 
page and hopefully it will work. Otherwise tell him again.


**Downloading files from uppmax to your computer**

To download files from uppmax to your local computer use scp
As an example when I download all the bamfiles from uppmax to my computer I write::
 scp johanr@milou.uppmax.uu.se:/proj/g2014046/webexport/files/RNAseqWorkshop/download/RNAseq/mappingsToRAB11FIP5/*.ba* .

 
**Creating your own download export page**

`Link to uppmax webexport guide 
<http://www.uppmax.uu.se/webexport-guide>`_

**RPKM multi mapping problem using multo**

`Link to multo 
<http://sandberg.cmb.ki.se/multo/>`_


**Answer regarding STAR**

`Link to forum where this is discussed and a quote from the creator of STAR  
<https://www.biostars.org/p/93883/>`_

qoute from Dobin
"the "Overhang" in these parameters has different meanings - bad naming choice, 
unfortunately.

The --sjdbOverhang is used only at the genome generation step, and tells STAR 
how many bases to concatenate from donor and acceptor sides of the junctions. 
If you have 100b reads, the ideal value of --sjdbOverhang is 99, which allows 
the 100b read to map 99b on one side, 1b on the other side. One can think of 
--sjdbOverhang as the maximum possible overhang for your reads.

On the other hand, --alignSJDBoverhangMin is used at the mapping step to define 
the minimum allowed overhang over splice junctions. For example, the default 
value of 3 would prohibit overhangs of 1b or 2b."





Lecture slides
==============
**Annotation of eukaryote genomes and transcriptomes, Henrik Lantz**

`Henrik Lantz 
<https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/lectures/Henrik_Lantz_Annotation.pdf>`_

**Transcriptome and isoform reconstruction with long reads, Ignas Bunikis**

Not present yet

**Proteogenomics, Janne Lehtiö**

Not present yet

**RNAseq QC, Johan Reimegård**

`Johan Reimegård 
<https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/lectures/Johan_Reimegard_RNAQC.pdf>`_


**Differential Expression, Mikael Huss**

`Mikael Huss 
<https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/lectures/Mikael_Huss_DiffExp.pdf>`_


**Differential expression in single-cell RNAseq, Åsa Björklund **

`Asa Björklund 
<https://export.uppmax.uu.se/g2014046/files/RNAseqWorkshop/download/RNAseq/lectures/Asa_Bjorklund_Single_Cell.pdf>`_



