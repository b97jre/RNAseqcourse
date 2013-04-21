==========================================
Assembling reads with Velvet
==========================================

Velvet
======
This assembler takes your reads and turns them into contigs. It consists of two
steps. In the first step, ``velveth``, the graph is created. Afterwards the
graph is traversed and contigs are created with ``velvetg``.

velveth
=======
Create the graph data structure with ``velveth``. We want to create assemblies
over multiple kmers. Again like we did with ``sickle``, first create a
directory with symbolic links to the pairs that you want to use::

    mkdir -p ~/glob/asm-workshop/velvet
    cd ~/glob/asm-workshop/velvet
    ln -s ../sickle/qtrim1.fastq pair1.fastq
    ln -s ../sickle/qtrim2.fastq pair2.fastq

Then the reads need to be interleaved for ``velveth``::

    shuffleSequences_fastq.pl pair1.fastq pair2.fastq pair.fastq

Run velveth over mutiple kmers::

    velveth out 21,58,6 -fastq -shortPaired pair.fastq

Check what directories have been created::

    ls

velvetg
=======
To get the actual contigs you will have to run ``velvetg`` over all of the created
graphs. You can vary options such expected coverage and the coverage cut-off if
you want, but we do not do that in this tutorial. We only choose not to do
scaffolding since we found that doesn't work that well::

    for dir in out_*; do velvetg $dir -scaffolding no; done

This runs ``velvetg`` sequentially on each directory. If you want to do it in
parallel instead running velvetg on each core, you can call ``velvetg`` using GNU
parallel ::

    echo out_* | sed 's/\ /\n/g' | parallel velvetg {} -scaffolding no

This changes the spaces from the echo output to newlines using ``sed``.  GNU
parallel then uses each line as an argument and places it at the location of
``{}``.

assemstats
==========
After the assembly one wants to look at the length distributions of the
resulting assemblies. You can use the ``assemstats`` script for that::

    assemstats 100 out_*/contigs.fa

I usually look at the distributions setting the cut-off at 100 and 1000 and put
the results in google docs. A useful program for copying the results from
assemstats is ``xclip``. This can copy the results to your clipboard so you can
directly past the results into google docs::

    assemstats 100 out_*/contigs.fa | xclip -selection clipboard


Note that ``xclip`` requires you to log in with X-forwarding enabled, i.e.
``ssh -X inod@kalkyl1.uppmax.uu.se``.
