==========================================
Merging reads with minimus2
==========================================

minimus2
========
Minimus2 is a script that uses programs from the Amos set of tools to merge
assemblies (http://sourceforge.net/apps/mediawiki/amos/index.php?title=AMOS).
Amos is quite a pain to install but luckily all tools are available again on
``/proj/b2010008/bin``.


I have created a script that basically follows the tutorial here:
http://ged.msu.edu/angus/metag-assembly-2011/velvet-multik.html. The script is
at ``/proj/b2010008/bin/merge-asm-minimus2.sh``. Please take a look at it and
try to understand what is going on::    

    less /proj/b2010008/bin/merge-asm-minimus2

It is basically a wrapper script that sets
an environment variable and calls the real script in my directory. Take a look
at that as well::

    less /bubo/home/h16/inod/glob/github/metassemble/scripts/assembly/merge-asm-minimus2.sh

And check the parameters to run it::

    merge-asm-minimus2 -h

Merge assemblies
=======================
Now we are going to try to merge assemblies from the previous velvet tutorial.
Go to the velvet directory and use ``merge-asm-minimus2``::
    
    cd ~/glob/asm-workshop/velvet
    merge-asm-minimus2 minimus2 out_*/contigs.fa

This gives you a file ``minimus2/all-merged.fasta`` that includes both the
merged and unmerged contigs. Do the same for the Ray assemblies::

    cd ~/glob/asm-workshop/ray
    merge-asm-minimus2 minimus2 out_*/Contigs.fasta

assemstats
==========
Compare all the results with assemstats::
    
    cd ~/glob/asm-workshop
    find | grep -e 'contigs.fa\|Contigs.fasta\|all-merged.fasta' | sort | xargs assemstats 100

``find`` finds all files from the current directory. Then grep uses a regular
expression (http://en.wikipedia.org/wiki/Regular_expression) to get the velvet
assemblies, Ray assemblies and minimus2 assemblies, sort them by name and give
the result as an argument to ``assemstats`` using ``xargs``.
