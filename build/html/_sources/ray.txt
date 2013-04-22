==========================================
Assembling reads with Ray
==========================================

Ray
===
Another assembler that takes your reads and turns them into contigs. It uses
the Message Passing Interface
(http://en.wikipedia.org/wiki/Message_Passing_Interface) to run over multiple
cores and nodes. Again ``Ray`` is located in ``/proj/b2010008/bin``.

First one has to load the ``openmpi`` module on uppmax::

    module load openmpi

Check if Ray is working with::

    mpiexec Ray --help

Ray goes a lot slower on one machine than velvet since it was optimized for
running on multiple. In this case the dataset is quite small so we will just be
using one machine anyway. Create the directory structure as we did before::
    
    mkdir -p ~/glob/asm-workshop/ray
    cd ~/glob/asm-workshop/ray
    ln -s ../velvet/pair.fastq pair.fastq

Then run Ray over multiple assemblies like we did with velvet. I recommend to
run this on an interactive node::

   for k in {45..57..6}; do mpiexec Ray -k $k -i pair.fastq -o out_$k; done

Note that using ``parallel`` here would not be benificial since you are using
multiple cores already with ``Ray``.

assemstats
==========
Look at the results with assemstats and compare them to the velvet output::

    cd ~/glob/asm-workshop/
    assemstats 100 ~/glob/asm-workshop/ray/out_*/Contigs.fasta ~/glob/asm-workshop/velvet/out_*/contigs.fa

