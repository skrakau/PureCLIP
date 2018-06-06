.. ` < >`_

Incorporating CL-motif scores
====================================

In order to address the crosslinking sequence bias, PureCLIP can incorporate information about motifs that are known to be preferentially crosslinked, also called CL-motifs (see also `Haberman et al., 2017 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1130-x>`_). 
For each CL-motif PureCLIP learns the influence on the crosslinking probability.
This can be particular useful for proteins binding to sequence motifs distinct from such CL-motifs.  

Compute CL-motif scores
-----------------------

To incorporate CL-motifs into the model of PureCLIP, first we need to compute position-wise CL-motif scores, indicating the positions CL-affinity.
You can use the provided precompiled list of `common_CL-motifs <https://github.com/skrakau/PureCLIP_data/blob/master/common_CL-motifs/>`_:   

.. code:: bash

    wget -O motifs.txt https://raw.githubusercontent.com/skrakau/PureCLIP_data/master/common_CL-motifs/dreme.w10.k4.txt
    wget -O motifs.xml https://raw.githubusercontent.com/skrakau/PureCLIP_data/master/common_CL-motifs/dreme.w10.k4.xml

If you want to compute the CL-motifs specific to the used eCLIP experiment, please have a look :ref:`here <customCLmotifs-label>`. 
Assume we have given a set of CL-motifs, we use `FIMO <http://meme-suite.org/doc/fimo.html>`_ (`Grant et. al, 2011 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3065696/>`_) to compute motif occurrences associated with a score within the reference sequence.
The following provided script (distributed with PureCLIP) first retrieves reference regions covered by the target experiment, then runs FIMO to compute position-wise CL-motif match scores within such regions and chooses for each position the motif with the highest score:

.. code:: bash

    export BEDTOOLS=/path/to/bedtools        # if not specified, PATH is searched
    export FIMO=/path/to/fimo                # if not specified, PATH is searched
    export WINEXTRACT=/path/to/winextract    # built together with PureCLIP
    compute_CLmotif_scores.sh ref.hg19.fa aligned.prepro.R2.bam motifs.xml motifs.txt fimo_clmotif_occurences.bed 

The resulting CL-motif scores are written to ``fimo_clmotif_occurences.bed``.       



PureCLIP
--------

The computed scores are then handed over to PureCLIP together with the parameter ``-nim 4``, indicating that scores with associated motif IDs 1-4 will be used (default: only scores with motif ID 1 are used). 

    .. code:: bash

        pureclip -i aligned.prepro.R2.bam -bai aligned.prepro.R2.bam.bai -g ref.hg19.fa -o PureCLIP.crosslink_sites.cov_CLmotifs.bed -nt 10 -iv 'chr1;chr2;chr3;' -nim 4 -fis fimo_clmotif_occurences.bed




