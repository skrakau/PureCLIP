
Incorporating input control data
====================================

.. container:: toggle

    .. container:: header

        **Generate sample files for minimal example: Show/Hide**

    |
   
    Beside the target PUM2 eCLIP data, we download the corresponding preprocessed `input data from ENCODE <https://www.encodeproject.org/experiments/ENCSR439GXW/>`_, and again filter the paired-end data to keep only R2:   

    .. code:: bash

       wget -O input.aligned.prepro.bam https://www.encodeproject.org/files/ENCFF043ERY/@@download/ENCFF043ERY.bam
       samtools view -hb -f 130 input.aligned.prepro.bam -o input.aligned.prepro.R2.bam
       samtools index input.aligned.prepro.R2.bam  

    
|

PureCLIP
--------

To run PureCLIP with input control data, additionally hand over the (preprocessed) BAM file from the input experiment with ``-ibam`` and the associated BAI file with ``-ibai``:

.. code:: bash

    pureclip -i aligned.prepro.R2.bam -bai aligned.prepro.R2.bam.bai -g ref.hg19.fa -iv 'chr1;chr2;chr3;' -nt 10 -o PureCLIP.crosslink_sites.cov_inputSignal.bed -ibam input.aligned.prepro.R2.bam -ibai input.aligned.prepro.R2.bam.bai



