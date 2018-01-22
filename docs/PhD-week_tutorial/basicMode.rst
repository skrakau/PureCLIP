.. PureCLIP documentation master file, created by
   sphinx-quickstart on Fri Jun 23 12:15:25 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Basic mode
====================================

.. container:: toggle

    .. container:: header

        **Get intermediate files: Show/Hide**

    |

    In case something went wrong, you can copy the intermediate files as follows:

    .. code:: bash

       cp ... 


    
|

PureCLIP
--------

To run PureCLIP in basic mode, it requires BAM and BAI files, the reference genome and a specified output file: 

.. code:: bash

    pureclip -i aligned.prepro.R2.bam -bai aligned.prepro.R2.bam.bai -g ref.hg19.fa -iv 'chr21;' -nt 8 -o PureCLIP.crosslink_sites.bed


With ``-iv`` the chromosomes (or transcripts) can be specified that will be used to learn the parameters of PureCLIPs HMM.
This reduces the memory consumption and runtime.
Usually, learning on a small subset of the chromosomes, e.g. Chr1-3, does not impair the results noticeable.
However, in the case of very sparse data this can be adjusted.
With ``-nt`` the number threads for parallelization can be specified. 


TODO
- bedtools -> expand regions
