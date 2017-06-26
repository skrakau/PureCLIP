
PureCLIPs output
====================================

Crosslink sites
---------------

The main output of PureCLIP is a BED6 file containing individual crosslink sites associated with a score:

1. **chr:** Name of the chromosome or scaffold.
2. **start:** Position of crosslink site.
3. **end:** Position behind crosslink site (start+1).
4. **state:** '3'
5. **score:** log posterior probability ratio of the first and second likely state.
6. **strand:** + or -


Binding regions
---------------

Optionally, if an output file for binding regions is specified with ``-or``, individual crosslink sites with a distance <= d (specified with ``-dm``, default 8 bp) are merged and given out in a separate BED6 file:
 
1. **chr:** Name of the chromosome or scaffold.
2. **start:** Start position, position of first crosslink site.
3. **end:** End position, position behind last crosslink site.
4. **indiv. scores:** 'score1;score2;score3;' 
5. **score:** Sum of log posterior probability ratio scores.
6. **strand:** + or -



