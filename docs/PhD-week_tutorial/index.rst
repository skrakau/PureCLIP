.. PureCLIP documentation master file, created by
   sphinx-quickstart on Fri Jun 23 12:15:25 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PhD Week '18 - Tutorial
====================================

In the following tutorial we will describe for each PureCLIP mode how to run a minimal example.


.. Hint::
    Since PureCLIP includes information from the transcriptomic neighbourhood, it is important to use a suitable reference sequence when mapping the reads. 
    Thus, when analysing CLIP data from proteins known to bind near exon-exon junctions on mRNAs, reads should be mapped directly against transcripts (e.g. as done by `Haberman et al., 2017 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1130-x>`_).

.. toctree::
   :hidden:
   :maxdepth: 2
   
   preprocessing
   basicMode
   incorporateInputAndCLmotifs
   sshmm


