.. PureCLIP documentation master file, created by
   sphinx-quickstart on Fri Jun 23 12:15:25 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Basic mode
====================================

.. rst-class:: html-toggle


.. container:: toggle

    .. container:: header

        **Show/Hide: Generate sample files for minimal example**
 
    As a first example you can download preprocessed `data from ENCODE <https://www.encodeproject.org/experiments/ENCSR661ICQ/>`_, and filter the paired-end data to keep only R2:

    .. code-block:: xml

       wget -O aligned.prepro.bam https://www.encodeproject.org/files/ENCFF280ONP/@@download/ENCFF280ONP.bam
       samtools view -hb -f 130 aligned.prepro.bam -o aligned.prepro.R2.bam
       samtools index aligned.prepro.R2.bam    

    Additionally, we need the corresponding reference genome:

    .. code-block:: xml

       wget -O ref.hg19.fa.gz https://www.encodeproject.org/files/female.hg19/@@download/female.hg19.fasta.gz 
       gunzip ref.hg19.fa.gz
    
|

PureCLIP
--------

To run PureCLIP in basic mode, it requires BAM and BAI files, the reference genome and a specified output file: 

.. code-block:: xml

    pureclip -i aligned.prepro.R2.bam -bai aligned.prepro.R2.bam.bai -g ref.hg19.fa -iv 'chr1;chr2;chr3;' -nt 10 -o PureCLIP.crosslink_sites.bed


With ``-iv`` the chromosomes (or transcripts) can be specified that will be used to learn the parameters of PureCLIPs HMM.
This reduces the memory consumption and runtime.
Usually, learning on a small subset of the chromosomes, e.g. Chr1-3, does not impair the results noticeable.
However, in the case of very sparse data this can be adjusted.
With ``-nt`` the number threads for parallelization can be specified. 






