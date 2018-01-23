.. PureCLIP documentation master file, created by
   sphinx-quickstart on Fri Jun 23 12:15:25 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree`directive.
.. ` < >`_

Preprocessing iCLIP/eCLIP-seq reads
====================================

This tutorial will walk you through an example how to preprocess your iCLIP/eCLIP data before analysing it.


RBFOX2 eCLIP data 
-----------------
 
We use `RBFOX2 eCLIP sequencing data <https://www.encodeproject.org/experiments/ENCSR756CKJ/>`_ from ENCODE (`Van Nostrand et. al, 2016 <https://www.ncbi.nlm.nih.gov/pubmed/27018577>`_).
However, since we have limited time and memory resources within this course, we will work only with a small subset of the data mapping to chr1, chr2 and chr21. 
You can find this data in ``~/protein-RNA-interactions/RBFOX2_data/rep1/`` and ``~/protein-RNA-interactions/RBFOX2_data/rep2/``.

.. code:: bash

    ls  ~/protein-RNA-interactions/RBFOX2_data/rep1/


The human reference sequences (containing chr1, chr2, chr21) as well as annotations are located in ``~/protein-RNA-interactions/hg19_data_data/``. 

.. code:: bash

   ls -lh ~/protein-RNA-interactions/hg19_data_data/ 


Let's start with changing to the folder 

.. code:: bash

    cd  ~/protein-RNA-interactions/RBFOX2_data/


Please note, that the following steps need to be done for both replicates, while only being described for rep1.


.. container:: toggle

    .. container:: header

        **Additional task: Show/Hide**

    |

    .. Hint::

        To preprocess the corresponding input control data, apply the same steps to the reads located in ``~/protein-RNA-interactions/input_control_data/``.
        Since only one input control replicate is given, the pooling step at the end is omitted. 



    
|



Adaptor trimming
----------------

Possible adapter contaminations at 3' ends can be removed using the tool `cutadapt <https://github.com/marcelm/cutadapt>`_ `(Martin, 2011) <http://journal.embnet.org/index.php/embnetjournal/article/view/200>`_ (as done by the YEO lab in the ENCODE `processing pipeline <https://www.encodeproject.org/documents/dde0b669-0909-4f8b-946d-3cb9f35a6c52/@@download/attachment/eCLIP_analysisSOP_v1.P.pdf>`_), while discarding reads shorter than 18bp: 

.. code:: bash

    cutadapt -j 6 --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT  -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA  -A CTTGTAGATCGGAAG  -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG  -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT  -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o rep1/reads.R1.trimmed.fastq.gz -p rep1/reads.R2.trimmed.fastq.gz rep1/reads.R1.fastq.gz rep1/reads.R2.fastq.gz > rep1/cutadapt.log

For eCLIP data, it was suggested (`Van Nostrand et. al, 2016 <https://www.ncbi.nlm.nih.gov/pubmed/27018577>`_) to apply two rounds of adapter trimming to correct for possible double ligations events.

.. code:: bash

    cutadapt -j 6 --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o rep1/reads.R1.trimmed2.fastq.gz -p rep1/reads.R2.trimmed2.fastq.gz rep1/reads.R1.trimmed.fastq.gz rep1/reads.R2.trimmed.fastq.gz > rep1/cutadapt.2.log


Preparing read IDs for UMI-tools
--------------------------------

We remove PCR duplicates based on the mapping position and random barcodes using UMI-tools, which requires the read ID in the format ``@HISEQ:87:00000000_BARCODE read1``.
However, in the current format the barcode is located in front of the actual read ID: 

.. code:: bash

    zless rep1/reads.R1.trimmed2.fastq.gz


Therefore we change the location of the barcode within the read ID prior the mapping, using the following commands:

.. code:: bash

    gunzip -c rep1/reads.R1.trimmed2.fastq.gz > rep1/reads.R1.trimmed2.fastq
    gunzip -c rep1/reads.R2.trimmed2.fastq.gz > rep1/reads.R2.trimmed2.fastq
    awk 'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (10+3), 500) "_" substr($1, 2, 10) " " $2 }; substr($1, 1, 1) != "@" {print}; ' rep1/reads.R1.trimmed2.fastq  | gzip > rep1/reads.R1.trimmed2.bc.fastq.gz
    awk 'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (10+3), 500) "_" substr($1, 2, 10) " " $2 }; substr($1, 1, 1) != "@" {print}; ' rep1/reads.R2.trimmed2.fastq  | gzip > rep1/reads.R2.trimmed2.bc.fastq.gz

where the used barcode length is 10.
TODO use easier command!!!!


Read mapping with STAR
----------------------

CLIP-seq reads can be mapped with the RNA-seq read aligner `STAR <https://github.com/alexdobin/STAR>`_ (`Dobin et. al, 2013 <https://www.ncbi.nlm.nih.gov/pubmed/23104886>`_).
It allows to include genome annotations in order to enable the alignment against spliced transcripts.
First, we need a genome index, created based on the reference sequences and the annotation file.
However, the preparation of this index requires > 8 GB of memory. 
You find an already created index in ``~/protein-RNA-interactions/hg19_data_data/genome_index/``.


.. Note::

   In general you can prepare your own genome index as follows

   .. code:: bash

       STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genome_index/ --genomeFastaFiles ref.fa --sjdbGTFfile annotation.gtf --sjdbOverhang 49

Next, we map the reads (R1 and R2) against the indexed genome:

.. code:: bash

    mkdir -p rep1/STAR
    STAR --outSAMtype BAM SortedByCoordinate --runThreadN 6 --genomeDir ~/protein-RNA-interactions/hg19_data_data/genome_index/ --readFilesIn rep1/reads.R1.trimmed2.bc.fastq.gz rep1/reads.R2.trimmed2.bc.fastq.gz --readFilesCommand  zcat --outFilterMultimapNmax 1 --scoreDelOpen -1 --outFileNamePrefix rep1/STAR/ --alignEndsType EndToEnd 

The parameter ``--outFilterMultimapNmax 1`` ensures only uniquely mapping reads will be reported.
Furthermore, it is important to use the ``--alignEndsType EndToEnd`` setting, to ensure the mapping of the whole read.
The aligned reads will be written then to STAR/Aligned.sortedByCoord.out.bam .

.. Note::

   Due to time and memory constraints within this course and since we prefiltered already the FASTQ files, we map the reads here only against the correspdonding subset of the genome, i.e. chr1, chr2, and chr21.
   In general it is recommended to use an assembly containing scaffolds as reference.
   This enables us to filter out reads that map both against a main chromosome and against a scaffold (e.g. ribosomal RNA).



Filtering
---------

Then we filter the aligned reads with `samtools <http://www.htslib.org/doc/samtools.html>`_  to obtain only reads that are mapped in proper pairs (``-f 2``) (a detailed explanation of available flags you can find `here <https://broadinstitute.github.io/picard/explain-flags.html>`_):
 
.. code:: bash

    samtools view -hb -f 2 rep1/STAR/Aligned.sortedByCoord.out.bam -o rep1/STAR/Aligned.f.bam 
    
and create an index, which is required for the next step

.. code:: bash

    samtools index rep1/STAR/Aligned.f.bam  


PCR duplicate removal using UMI-tools
-------------------------------------

For truncation based CLIP-seq data it is crucial to remove PCR duplicates to allow for an accurate crosslink site detection.
We use the `UMI-tools <https://github.com/CGATOxford/UMI-tools>`_ (`Smith et. al, 2017 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/>`_), which is able to handle errors within barcode sequences.

.. code:: bash

    umi_tools dedup -I rep1/STAR/Aligned.f.bam --paired -S rep1/STAR/Aligned.f.duplRm.bam > rep1/STAR/umi_tools.log
 

Pooling and R2 retrieval
------------------------

Finally, we merge the preprocessed alignments of the individual replicates:

.. code:: bash

    samtools merge -f Aligned.f.duplRm.pooled.bam rep1/STAR/Aligned.f.duplRm.bam rep2/STAR/Aligned.f.duplRm.bam

and filter for R2, to keep only reads containing information about potential truncation events (for iCLIP data this would be R1):

.. code:: bash

    samtools view -hb -f 130 Aligned.f.duplRm.pooled.bam -o Aligned.f.duplRm.pooled.R2.bam
    samtools index Aligned.f.duplRm.pooled.R2.bam   


Quality control
---------------

It's always a good idea to assess the quality of the data prior to the actual analysis.
For this we use `fastqc <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_:

.. code:: bash

    mkdir fastqc
    fastqc -o fastqc/ aligned.f.duplRm.pooled.R2.bam



