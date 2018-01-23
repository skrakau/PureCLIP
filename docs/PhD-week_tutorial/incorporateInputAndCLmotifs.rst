.. ` < >`_

PureCLIP: incorporating input control data and CL-motifs
=========================================================



Input control signal
-----------------------

In case you did not already preprocess the input control data, you can get the intermediate BAM and BAI files as follows:
   
.. code:: bash

    cp ../intermediate_results/input_control_data/Aligned.f.duplRm.R2.bam  ../input_control_data/
    samtools index ./input_control_data/Aligned.f.duplRm.R2.bam


CL-motifs 
------------------------

In order to address the crosslinking sequence bias, PureCLIP can incorporate information about motifs that are known to be preferentially crosslinked, also called CL-motifs (see also `Haberman et al., 2017 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1130-x>`_). 
For each CL-motif PureCLIP learns the influence on the crosslinking probability.
This can be particular useful for proteins binding to sequence motifs distinct from such CL-motifs, e.g. RBFOX2. 

We need to know the positions of CL-motif occurrences (more details described in `incorporateCLmotifs.html <http://pureclip.readthedocs.io/en/latest/PureCLIPTutorial/incorporateCLmotifs.html>`_).
Here we use a given BED file ``~/protein-RNA-interactions/hg19_data/CLmotif_occurences.chr1_2_21.bed`` containing already computed occurrences of a set of four common CL-motifs together with a score.





PureCLIP
--------

To run PureCLIP with input control data, additionally hand over the (preprocessed) BAM file from the input experiment with ``-ibam`` and the associated BAI file with ``-ibai``.
The computed CL-motif occurrences are then handed over with ``-fis`` together with the parameter ``-nim 4``, indicating that scores with associated motif IDs 1-4 will be used (default: only scores with motif ID 1 are used). 

.. code:: bash

    pureclip -i Aligned.f.duplRm.pooled.R2.bam -bai Aligned.f.duplRm.pooled.R2.bam.bai -g ../hg19_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.chr1_2_21.fa -iv 'chr21;' -bdw 20 -nt 8 -ibam ../input_control_data/Aligned.f.duplRm.R2.bam -ibai ../input_control_data/Aligned.f.duplRm.R2.bam.bai -nim 4 -fis ../hg19_data/CLmotif_occurences.chr1_2_21.bed -o crosslinkSites.input_CLmotifs.bed -or bindingRegions.input_CLmotifs.bed > pureclip.input_CLmotifs.log





