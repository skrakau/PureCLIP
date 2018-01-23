.. PureCLIP documentation master file, created by
   sphinx-quickstart on Fri Jun 23 12:15:25 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PureCLIP: basic mode
====================================

.. container:: toggle

    .. container:: header

        **Get intermediate files: Show/Hide**

    |

    In case something went wrong during the preprocessing, you can obtain the intermediate files as follows:

    .. code:: bash

       cp ~/protein-RNA-interactions/intermediate_results/RBFOX2_data/Aligned.f.duplRm.pooled.R2.bam .


    
|

PureCLIP
--------

To run PureCLIP in its basic mode, i.e. without incorporating external data as covariates, it requires BAM and BAI files, the reference genome and specified output files: 

.. code:: bash

    mkdir PureCLIP_results
    pureclip -i Aligned.f.duplRm.pooled.R2.bam -bai Aligned.f.duplRm.pooled.R2.bam.bai -g ~/protein-RNA-interactions/hg19_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr1_2_21.fa -iv 'chr21;' -bdw 20 -nt 8 -o PureCLIP_results/crosslinkSites.basic.bed -or PureCLIP_results/bindingRegions.basic.bed > PureCLIP_results/pureclip.basic.log


With ``-iv`` the chromosomes (or transcripts) can be specified that will be used to learn the parameters of PureCLIPs HMM.
This reduces the memory consumption and runtime.
Usually, learning on a small subset of the chromosomes, e.g. Chr1-3, does not impair the results noticeable.
However, in the case of very sparse data this can be adjusted.
With ``-nt`` the number threads for parallelization can be specified. 
The parameter ``-bdw`` is the bandwidth used for smoothing to read start counts (default: ``-bdw 50``).   


