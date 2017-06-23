.. PureCLIP documentation master file, created by
   sphinx-quickstart on Fri Jun 23 12:15:25 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree`directive.
.. ` < >`_

Preprocessing iCLIP/eCLIP-seq reads
====================================

This tutorial will walk you through an example how to preprocess your iCLIP/eCLIP data before analysing it.
This only serves as a guide and you can of course as well use your own CLIP-seq preprocessing workflow.

Requirements for this tutorial
-------------------------------

The tutorial requires the following tools (and was tested for the specified versions):

* `cutadapt <https://cutadapt.readthedocs.io/en/stable/installation.html>`_ (v1.12)
* `samtools <https://github.com/samtools/samtools>`_ (v0.1.19-44428cd)
* `STAR  <https://github.com/alexdobin/STAR>`_ (v2.5.1b)
* `UMI <https://github.com/CGATOxford/UMI-tools>`_ (v0.4.3)
* `fastqc <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ (v0.11.5)


Data retrieval
--------------

Let's start with downloading the data.
We use the raw `PUM2 eCLIP sequencing data <https://www.encodeproject.org/experiments/ENCSR661ICQ/>`_ from ENCODE (`Van Nostrand et. al, 2016 <https://www.ncbi.nlm.nih.gov/pubmed/27018577>`_):

.. code:: bash

    mkdir rep1 rep2 
    wget -O rep1/reads.R1.fastq.gz https://www.encodeproject.org/files/ENCFF956TOZ/@@download/ENCFF956TOZ.fastq.gz
    wget -O rep1/reads.R2.fastq.gz https://www.encodeproject.org/files/ENCFF133DNM/@@download/ENCFF133DNM.fastq.gz
    wget -O rep2/reads.R1.fastq.gz https://www.encodeproject.org/files/ENCFF041KJT/@@download/ENCFF041KJT.fastq.gz
    wget -O rep2/reads.R2.fastq.gz https://www.encodeproject.org/files/ENCFF462SCV/@@download/ENCFF462SCV.fastq.gz

We download the human reference genome (primary assembly containing chromosomes and scaffolds) and the genome annotation (comprehensive set) from `GENCODE <http://www.gencodegenes.org/>`_:

.. code:: bash

    wget -O ref.GRCh38.fa.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz 
    gunzip ref.GRCh38.fa.gz

    wget -O annotation.v26.gtf.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.primary_assembly.annotation.gtf.gz
    gunzip annotation.v26.gtf.gz

Please note, that the following steps need to be done for both replicates, while only being described for rep1.


Adaptor trimming
----------------

Possible adapter contaminations at 3' ends can be removed using the tool `cutadapt <https://github.com/marcelm/cutadapt>`_ `(Martin, 2011) <http://journal.embnet.org/index.php/embnetjournal/article/view/200>`_ (as done by the YEO lab in the ENCODE `processing pipeline <https://www.encodeproject.org/documents/dde0b669-0909-4f8b-946d-3cb9f35a6c52/@@download/attachment/eCLIP_analysisSOP_v1.P.pdf>`_), while discarding reads shorter than 18bp: 

.. code:: bash

    cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT  -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA  -A CTTGTAGATCGGAAG  -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG  -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT  -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o rep1/reads.R1.trimmed.fastq.gz -p rep1/reads.R2.trimmed.fastq.gz rep1/reads.R1.fastq.gz rep1/reads.R2.fastq.gz

For eCLIP data, it was suggested (`Van Nostrand et. al, 2016 <https://www.ncbi.nlm.nih.gov/pubmed/27018577>`_) to apply two rounds of adapter trimming to correct for possible double ligations events.

.. code:: bash

    cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o rep1/reads.R1.trimmed2.fastq.gz -p rep1/reads.R2.trimmed2.fastq.gz rep1/reads.R1.trimmed.fastq.gz rep1/reads.R2.trimmed.fastq.gz


Preparing read IDs for UMI
--------------------------

We remove PCR duplicates based on the mapping position and random barcodes using UMI, which requires the read ID in the format ``_@HISEQ:87:00000000_BARCODE read1_``.
Therefore we append the barcode to the read ID prior the mapping, using the following command:

.. code:: bash

    zcat rep1/reads.R1.trimmed2.fastq.gz > rep1/reads.R1.trimmed2.fastq
    zcat rep1/reads.R2.trimmed2.fastq.gz > rep1/reads.R2.trimmed2.fastq
    awk -v l=10 'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }; substr($1, 1, 1) != "@" {print}; ' rep1/reads.R1.trimmed2.fastq  | gzip > rep1/reads.R1.trimmed2.bc.fastq.gz
    awk -v l=10 'BEGIN{OFS=FS=" "} substr($1, 1, 1) == "@" {print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }; substr($1, 1, 1) != "@" {print}; ' rep1/reads.R2.trimmed2.fastq  | gzip > rep1/reads.R2.trimmed2.bc.fastq.gz

where ``l=10`` denotes the used barcode length.


Read mapping with STAR
----------------------

CLIP-seq reads can be mapped with the RNA-seq read aligner `STAR <https://github.com/alexdobin/STAR>`_ (`Dobin et. al, 2013 <https://www.ncbi.nlm.nih.gov/pubmed/23104886>`_).
It allows to include genome annotations in order to enable the alignment against spliced transcripts.
First, we need to prepare the genome index, using the annotation file:

.. code:: bash

    mkdir genome_index
    STAR --runThreadN 10 --runMode genomeGenerate --genomeDir genome_index/ --genomeFastaFiles ref.GRCh38.fa --sjdbGTFfile annotation.v26.gtf --sjdbOverhang 49

Next, we map the reads (R1 and R2) against the indexed genome:

.. code:: bash

    mkdir -p rep1/STAR
    STAR --outSAMtype BAM SortedByCoordinate --runThreadN 10 --genomeDir genome_index/ --readFilesIn rep1/reads.R1.trimmed2.bc.fastq.gz rep1/reads.R2.trimmed2.bc.fastq.gz --readFilesCommand  zcat --outFilterType BySJout --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --scoreDelOpen -1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix rep1/STAR/ --alignEndsType EndToEnd 

The parameter ``--outFilterMultimapNmax 1`` ensures only uniquely mapping reads will be reported.
Since we used the primary assembly containing scaffolds as reference, this enables us to filter out reads that map both against a main chromosome and against a scaffold (e.g. ribosomal RNA).
Furthermore, it is important to use the ``--alignEndsType EndToEnd`` setting, to ensure the mapping of the whole read.
The aligned reads will be written then to STAR/Aligned.sortedByCoord.out.bam .


Filtering
---------

We filter the aligned reads to obtain only reads mapping against the main chromosomes:
 
.. code:: bash

    samtools index rep1/STAR/Aligned.sortedByCoord.out.bam
    samtools view -hb -f 2 rep1/STAR/Aligned.sortedByCoord.out.bam -o rep1/aligned.f.bam chr1:1 chr2:1 chr3:1 chr4:1 chr5:1 chr6:1 chr7:1 chr8:1 chr9:1 chr10:1 chr11:1 chr12:1 chr13:1 chr14:1 chr15:1 chr16:1 chr17:1 chr18:1 chr19:1 chr20:1 chr21:1 chr22:1 chrX:1 chrY:1
    samtools index rep1/aligned.f.bam


PCR duplicate removal using UMI
--------------------------------

For truncation based CLIP-seq data it is crucial to remove PCR duplicates to allow for an accurate crosslink site detection.
We use the tool `UMI <https://github.com/CGATOxford/UMI-tools>`_ (`Smith et. al, 2017 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/>`_), which is able to handle errors within barcode sequences.

.. code:: bash

    umi_tools dedup -I rep1/aligned.f.bam --paired -S rep1/aligned.f.duplRm.bam
 

Pooling and R2 retrieval
------------------------

Finally, we merge the preprocessed alignments of the individual replicates:

.. code:: bash

    samtools merge -f aligned.f.duplRm.pooled.bam rep1/aligned.f.duplRm.bam rep2/aligned.f.duplRm.bam

and filter for R2, to keep only reads containing information about potential truncation events (for iCLIP data this would be R1):

.. code:: bash

    samtools view -hb -f 130 aligned.f.duplRm.pooled.bam -o aligned.f.duplRm.pooled.R2.bam
    samtools index aligned.f.duplRm.pooled.R2.bam


Quality control
---------------

It's always a good idea to assess the quality of the data prior to the actual analysis.
For this we use `fastqc <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_:

.. code:: bash

    mkdir fastqc
    fastqc -o fastqc/ aligned.f.duplRm.pooled.R2.bam

