# column4
## Initial condition
Magma column at thermodynamic equilibrium with
**total water at 2 wt%, total CO<sub>2</sub> at 0.5 wt%**, resulting in a
volatile-saturated magma. Excess fluids are removed.
**Silicate melt composition** is taken from [Métrich et al., 2010](https://doi.org/10.1093/petrology/egp083), sample
ST531 4b

| SiO2  | TiO2  | Al2O3 | Fe2O3 | FeO   | MnO   | MgO   | CaO   | Na2O | K2O  |
|-------|-------|-------|-------|-------|-------|-------|-------|------|------|
| .4741 | .0087 | .1496 | .0164 | .0532 | .0009 | .0394 | .1556 | .022 | .016 |

```
oxidesList = [.02, .005, .4741, .0087, .1496, .0164, .0532, .0009, .0394, .1556, .022, .016]
```

The system is isothermal, with temperature T = 1473 K. Pressure range for the column is 400
MPa to 0.1 MPa (atmospheric pressure), and pressure steps for
thermodyamic equilibrium calculations are at 10 MPa.

## Flushing
The flushing fluid phase is injected in discrete batches; every batch
is 1 wt% of fluid with respect to magma. The flushing fluid
composition is 0.94 mol CO2. At every
pressure step equilibrium is recalculated, the excess fluid is removed
and passed on to the next pressure step.

## Results
Outputs are matrix text files in [results](pyResults), one file for
each input batch as indicated in the file name. Columns are

| p | T | total water | total CO2 | dissolved water | dissolved CO2 | exsoved water | exsolved water mol | fluid |
|---|---|-------------|-----------|-----------------|---------------|---------------|--------------------|-------|
|   |   |             |           |                 |               |               |                    |       |

All data are expressed as weight fractions except where otherwise
noted. Dissolved quantities are with respect to total melt phase,
while exsolved quantities are with respect to total fluid phase (water
and CO2).
  
## Credits
Chiara P Montagna, INGV Pisa  
chiara.montagna@ingv.it
