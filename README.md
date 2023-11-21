# pyFlushing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10176157.svg)](https://doi.org/10.5281/zenodo.10176157)

A python main to simulate CO<sub>2</sub> flushing in a magma column through equilibrium steps, based on [solwcad](http://www.pi.ingv.it/progetti/eurovolc/#SOLWCAD) [(Papale et al. 2006)](https://doi.org/10.1016/j.chemgeo.2006.01.013).
Needs [`numpy`](https://numpy.org/).

First, `solwcad` needs to be wrapped in python. To do so, use
[`f2py`](https://numpy.org/doc/stable/f2py/):

```
python -m numpy.f2py -c -m solwcadsub solwcadsub.f
```

Two calculation modes are available:

* [flushing column](flushingColumn)
Initial condition is a magma column with given total volatile (water
and CO<sub>2</sub>) content at equilibrium, with no excess fluids. The
flushing fluid phase is injected in discrete batches. At every
pressure step equilibrium is recalculated, the excess fluid is removed
and passed on to the next pressure step, while the melt composition changes before the next flushing batch enters from below.

* [flushing equilibrium](flushingEquilibrium)
Initial condition is a magma column in equilibrium with a
fixed-composition fluid phase (isopleth), with no excess fluids. The
flushing happens as above.

Simulated cases are described in their respective folders
[stromboliLowCO2](flushingColumn/stromboliLowCO2) and
[stromboliHighCO2](flushingColumn/stromboliHighCO2).

Scripts for plotting results using [gnuplot](http://www.gnuplot.info/)
and [R](https://www.r-project.org/) are also provided.

### References
Papale, P., Moretti, R., & Barbato, D. (2006). The compositional dependence of the saturation surface of H2O+CO2 fluids in silicate melts. Chemical Geology, 229(1–3), 78–95. https://doi.org/10.1016/j.chemgeo.2006.01.013

## Credits
Chiara P Montagna, INGV Pisa 
chiara.montagna@ingv.it
