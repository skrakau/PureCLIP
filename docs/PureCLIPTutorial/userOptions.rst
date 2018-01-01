
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




