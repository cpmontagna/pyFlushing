# column1
## Initial condition
Magma column at thermodynamic equilibrium with
**total water at 2 wt%, total CO<sub>2</sub> at 0.5 wt%**, resulting in a
volatile-saturated magma. Excess fluids are removed.
**Silicate melt composition** is taken from Metrich et al., 2010; sample
ST531 4b

| SiO2  | TiO2  | Al2O3 | Fe2O3 | FeO   | MnO   | MgO   | CaO   | Na2O | K2O  |
|-------|-------|-------|-------|-------|-------|-------|-------|------|------|
| .4741 | .0087 | .1496 | .0164 | .0532 | .0009 | .0394 | .1556 | .022 | .016 |

```
oxidesList = [.02, .005, .4741, .0087, .1496, .0164, .0532, .0009, .0394, .1556, .022, .016]
```

## Flushing
The flushing fluid phase is injected in discrete batches. At every
pressure step equilibrium is recalculated, the excess fluid is removed
and passed on to the next pressure step.



Chiara P Montagna, INGV Pisa, 25/2/2021
chiara.montagna@ingv.it
