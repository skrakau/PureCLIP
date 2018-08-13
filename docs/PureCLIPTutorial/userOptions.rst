
General user options
====================================


For a full list of user options, please type:


.. code:: bash

    pureclip --help



.. Hint::
    By default, the parameters are set to values optimized for proteins binding to short defined binding regions, e.g. proteins binding to short specific motifs such as PUM2 and RBFOX2.
    With the ``-bc``  option this behaviour can be changed.
    Use ``-bc 1`` for proteins causing larger crosslink clusters with relatively lower read start counts (e.g. proteins binding to low complexity motifs).

.. Hint::
    In case of data suffering from a large fraction of non-coinciding read start counts, for example caused by constrained ends, decrease the initial binomial probability parameter of the 'crosslink' state (e.g. ``-b2p 0.03``). 
    Additionally use the option ``-antp`` to enable the computation of a suitable n threshold, which is applied when learning parameters linked to the 'crosslink' state.
    With this we avoid the parameter learning to be impaired by lower covered regions where no distinction between 'crosslink' and 'non-crosslink' states is possible.


Remarks
---------------

Floating point precision:

 - By default PureCLIP stores emission, state posterior and forward-backward probabilities as double floating-point numbers. However, for some data, in particular when including background control data that additionally contains artefacts, emission probabilities of outliers can become very small. In such cases, in order to allow the required computations, a higher floating point precision is required, i.e. long double precision, which can be applied with the parameter ``-ld``. Note that this comes with a higher memory consumption.


Training set to learn model parameters:

 - In order to reduce the memory consumption of PureCLIP, we learned the model parameters used in the PureCLIP paper only for a subset of chromosomes, i.e. ``-iv 'chr1;chr2;chr3;'``. When using PureCLIP in basic mode, i.e. without incorporating any covariates, for the evaluated data this does not significantly change the results. However, it should be noted that when incorporating input signal, PureCLIPs precision usually improves when learning on a larger set.


Gamma shape parameters when incorporating background control data:

 - By default, PureCLIP enforces the shape parameter of the 'non-enriched' gamma distribution to be less or equal than the shape parameter of the 'enriched' distribution. This constraint can be turned off using ``-fk`` (it could be observed to improve results for some datasets).


Mapping artefacts:

 - In general mapping artefacts should be handled during preprocessing. However, the parameter ``-mkn`` can be used to define the max. k/N ratio (#read start sites/fragment coverage in region) when learning truncation probabilities for the 'non-crosslink' and 'crosslink' states. A large number of extreme high ratios, e.g. 1, might originate from mapping artifacts and can disturb the parameter learning so that PureCLIP becomes insensitive for sites with lower ratios.

