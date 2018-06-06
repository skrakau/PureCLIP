.. _customCLmotifs-label:


Compute custom CL-motifs
------------------------

1.  Detect crosslink sites on the input data using the basic version of PureCLIP.


2.  Learn CL-motifs by running `DREME <http://meme-suite.org/doc/dreme.html>`_ (`Bailey et. al, 2011 <https://www.ncbi.nlm.nih.gov/pubmed/21543442>`_) with the parameters ``-norc -k 6 -m 4`` on 10-bp windows spanning the called input crosslink sites, while using 10-bp windows 20 bp upstream and downstream as the control.

  
