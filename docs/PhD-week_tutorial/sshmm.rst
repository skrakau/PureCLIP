.. 

ssHMM: extracting sequence-structure motifs
===========================================

.. container:: toggle

    .. container:: header

        **Get intermediate files: Show/Hide**

    |

    In case something went wrong, you can copy the intermediate files as follows:

    .. code:: bash

       cp ... 


    
|

Preprocessing
-------------

We need to extend PureCLIP regions +- 20 bp and get the top 3000 ... bla

.. code:: bash

    bedtools ?
    sort -g -k5,5 -r PureCLIP_results/regions.basic.bed > PureCLIP_results/regions.basic.soSc.bed
    head -n 3000 PureCLIP_results/regions.basic.soSc.bed > PureCLIP_results/regions.basic.top3000.bed

To use the binding regions for ssHMM, we copy these to ...
    
.. code:: bash

    cp PureCLIP_results/regions.basic.top3000.bed datasets/bed/RBFOX2_PureCLIP-basic_regions/positive_raw.bed


Blabla

.. code:: bash

    preprocess_dataset -e 20 datasets RBFOX2_PureCLIP-basic_regions 1 0.0

ssHMM
-------------

Blabla

.. code:: bash

    train_seqstructhmm datasets/fasta/RBFOX2_PureCLIP-basic_top3000_regions/positive.fasta datasets/shapes/RBFOX2_PureCLIP-basic_top3000_regions/positive.txt -o results/ -n 6 -b -j RBFOX2_PureCLIP-basic_top3000_regions_len6_b_random


