.. ` < >`_

Incorporating input control data and CL-motif scores
====================================================

In order to address the crosslinking sequence bias, PureCLIP can incorporate information about motifs that are known to be preferentially crosslinked, also called CL-motifs (see also `Haberman et al., 2017 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1130-x>`_). 
For each CL-motif PureCLIP learns the influence on the crosslinking probability.
This can be particular useful for proteins binding to sequence motifs distinct from such CL-motifs.  

Get CL-motif scores
-----------------------
  

You can get them here


PureCLIP
--------

The computed scores are then handed over to PureCLIP together with the parameter ``-nim 4``, indicating that scores with associated motif IDs 1-4 will be used (default: only scores with motif ID 1 are used). 

    .. code:: bash

        pureclip -i aligned.prepro.R2.bam -bai aligned.prepro.R2.bam.bai -g ref.hg19.fa -o PureCLIP.crosslink_sites.cov_CLmotifs.bed -nt 10 -iv 'chr1;chr2;chr3;' -nim 4 -fis fimo_clmotif_occurences.bed




