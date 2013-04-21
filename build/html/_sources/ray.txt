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
using one machine anyway. I recommend to run this on an interactive node.
Create the directory structure as we did before::
    
    mkdir -p ~/glob/asm-workshop/ray
    cd ~/glob/asm-workshop/ray
    ln -s ../velvet/pair.fastq pair.fastq

Then run Ray::

   mpiexec 
