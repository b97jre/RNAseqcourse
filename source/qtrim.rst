==========================================
Quality trimming Illumina paired-end reads
==========================================

Sickle
======
For quality trimming Illumina paired end reads we use the library sickle which
trims reads from 3' end to 5' end using a sliding window. If the mean quality
drops below a specified number remaining part of the read will be trimmed.


Sickle on Uppmax
================
Daniel has already installed sickle in the project dir on uppmax so hurray! The
absolute path is ``/proj/b2010008/bin/sickle``. ``/proj/b2010008/bin``` to your
``PATH`` variable if you haven't done so already::

    export PATH=$PATH:/proj/b2010008/bin

Note that this sets the variable only for your current shell and all
subprocesses spawned by it. If you want the directory added to your path for
always add the command to the file ``~/.bashrc`` which is executed everytime bash
is started.

(Optional) Install sickle by yourself
=====================================
Follow these steps only if you want to install sickle by yourself. We are not
doing this for the workshop today, but I included it just for future reference.
Installation procedures of research software are often follow the same pattern.
Download the code, *compile* it and copy the binary to a location in your
``PATH``. The code for sickle is on https://github.com/najoshi/sickle. I prefer
*compiling* my programs in ``~/src`` and then copying the resulting program to
my ``~/bin`` directory, which is in my ``PATH``. This should get you a long
way::

    mkdir -p ~/src

    # Go to the source directory and clone the sickle repository
    cd ~/src
    git clone https://github.com/najoshi/sickle
    cd sickle

    # Compile the program
    make

    # Create a bin directory
    mkdir -p ~/bin
    cp sickle ~/bin


Downloading a test set
======================
Today we'll be working on a small metagenomic data set from the anterior nares
(http://en.wikipedia.org/wiki/Anterior_nares).

.. image:: ../images/nostril.jpg


So get ready for your first smell of metagenomic assembly - O, yes I did. Run
all these commands in your shell::
    
    # Download the reads and extract them
    mkdir -p ~/glob/asm-workshop/data
    cd ~/glob/asm-workshop/data
    wget http://downloads.hmpdacc.org/data/Illumina/anterior_nares/SRS018585.tar.bz2
    tar -xjf SRS018585.tar.bz2

If successfull you should have the files
``~/asm-workshop/data/SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq``
and
``~/asm-workshop/data/SRS018585/SRS018585.denovo_duplicates_marked.trimmed.2.fastq``.


Running sickle on a paired end library
======================================
I like to create directories for specific parts I'm working on and creating
symbolic links (shortcuts in windows) to the files so there we go::
    
    mkdir -p ~/glob/asm-workshop/sickle
    cd ~/glob/asm-workshop/sickle
    ln -s ../data/SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq pair1.fastq
    ln -s ../data/SRS018585/SRS018585.denovo_duplicates_marked.trimmed.1.fastq pair2.fastq

The difference between a symbolic link a and a hard link can be found here:
http://stackoverflow.com/questions/185899/what-is-the-difference-between-a-symbolic-link-and-a-hard-link.
In this case I use symbolic links so I know what path the original reads have,
which can help one remember what those reads were after an intoxicating night
out. Now run sickle::

    # check if sickle is in your PATH
    which sickle
    # Run sickle
    sickle pe \
        -f pair1.fastq \
        -r pair2.fastq \
        -t sanger \
        -o qtrim1.fastq \
        -p qtrim2.fastq \
        -s qtrim.unpaired.fastq
    # Check what files have been generated
    ls
