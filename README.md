A small python library to simulate flushing of a magma column through equilibrium steps, based on [solwcad](http://www.pi.ingv.it/progetti/eurovolc/) [(Papale et al. 2006)](https://doi.org/10.1016/j.chemgeo.2006.01.013)

* Dependencies
Needs `numpy`.

* [flushing column](flushingColumn)
Initial condition is a magma column with given total volatile (water and CO<sub>2</sub>) content at equilibrium, with no excess fluids. The flushing fluid phase is injected in discrete batches. At every pressure step equilibrium is recalculated, the excess fluid is removed and passed on to the next pressure step, changing the composition of the melt before the next flushing batch is inserted.

* [flushing equilibrium](flushingEquilibrium)
Initial condition is a magma column in equilibrium with a fixed-composition fluid phase (isopleth), with no excess fluids. The flushing happens as above.

Chiara P Montagna, INGV Pisa, 25/2/2021
chiara.montagna@ingv.it
