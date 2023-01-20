# simulate flushing of a magma column through equilibrium steps
# initial condition: magma column with given total volatile content at equilibrium, with no excess fluids
# excess fluid phase is injected in discrete batches
# at every pressure step equilibrium is recalculated, the excess fluid is removed and passed on to the next pressure step, changing the composition of the melt before the next flushing batch is inserted
# equilibrium calculations performed with solwcad (Papale et al., 2006)

# INPUTS:
# silicate melt composition in terms of 10 major oxides
# total water and CO2 content of the initial magma column
# temperature (calculations are isothermal)
# pressure range of the magma column
# flushing fluid amount and composition

# OUTPUTS:
# a matrix text file for each flushing input step, with columns:
# temperature
# pressure
# total water
# total CO2
# dissolved water, weight fraction on total
# dissolved co2, weight fraction on total
# exsolved water, weight fraction on gas
# exsolved water, molar fraction on gas
# weight fraction of fluid
 
# needs numpy

# Chiara P Montagna, INGV Pisa, 25/2/2021
# chiara.montagna@ingv.it

import solwcadsub as sw        # import solwcad routine
import numpy as np
import os

# ============ INPUTS =========================================================

T = 1473.                      # temperature, K
initial_pressure = 1.e5        # magma column pressure range, Pa
final_pressure = 4.e8
delta_pressure = 1.e7

# magma composition: weight fractions of water, CO2 and ten major oxides (SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O)
# silicate melt composition: sample ST531 4b Metrich et al. 2010
# initial water at 2 wt%, initial CO2 at 0.5 wt%, volatile-saturated magma
oxidesList = [.02, .005, .4741, .0087, .1496, .0164, .0532, .0009, .0394, .1556, .022, .016]

# flushing phase
# injected in the system at steps delta_input_gas, from initial_input_gas to final_input_gas
initial_input_gas = 0.01
final_input_gas = 10.
delta_input_gas = 0.01 
co2_in = 0.97              # composition (weight fraction) of the fluid phase

# save results
saveDir = 'pyResults'    # directory where to save results
initialCondition = 'initialColumnPy.out' # file name for the initial condition
saveFile = 'pythonColumn'  # file names for successive flushing batches - total fluid flushed so far is appended

# ====================================================================

ox = np.array(oxidesList)

if not os.path.exists(saveDir):
    os.makedirs(saveDir)

# molar masses
M_h2o = 0.0018
M_co2 = 0.0044

# loops on input flushing gas and pressure
input_size = int((final_input_gas - initial_input_gas)/delta_input_gas) + 1
vector_size = int((final_pressure - initial_pressure)/delta_pressure) + 1

previous_dissolved_h2o = np.zeros(vector_size)
previous_dissolved_co2 = np.zeros(vector_size)

# initial condition: magma column in equilibrium with initial water and co2, with no excess fluid
with open(saveDir + '/' + initialCondition,'w') as fs:
    fs.write("weight fractions\n")
    fs.write("pressure \t temperature \t water \t CO2 \t dissolved h2o onliq \t dissolved co2 onliq \t exsolved h2o ongas \t exsolved h2o ongas mol \t gas_ontotal\n")

for i in range(0, vector_size):
    
    p = final_pressure - (delta_pressure * i)
    xy = sw.solwcad(p,T,ox,3)
    exsolved_h2o_ongas = xy[2]

    exsolved_h2o_ongas_molar = exsolved_h2o_ongas * M_co2 / (exsolved_h2o_ongas * (M_co2 - M_h2o) + M_h2o)

    gas_ontotal = (ox[0] + ox[1] - xy[0] - xy[1])/(1. - xy[0] - xy[1])

    with open(saveDir + '/' + initialCondition,'a') as fs:
        fs.write("%f \t %f \t %f \t %f \t % f \t %f \t %f \t %f \t %f \n" % (p, T, ox[0], ox[1], xy[0], xy[1], xy[2], exsolved_h2o_ongas_molar, gas_ontotal))

    total_exsolved_h2o = gas_ontotal * exsolved_h2o_ongas
    total_exsolved_co2 = gas_ontotal * (1. - exsolved_h2o_ongas)

    previous_dissolved_h2o[i] = xy[0] * (1. - gas_ontotal)
    previous_dissolved_co2[i] = xy[1] * (1. - gas_ontotal) 
      
    ox[0] = ox[0] - total_exsolved_h2o
    ox[1] = ox[1] - total_exsolved_co2

# flushing: at each pressure step, the excess fluid phase is added to
# dissolved volatiles in the next step

for j in range(0, input_size):
    
    input_gas = initial_input_gas + delta_input_gas * j
    with open(saveDir + '/' + saveFile + "{:.2f}".format(input_gas) + '.out','w') as fs:
    
        fs.write("weight fractions unless otherwise noted\n")
        fs.write("pressure \t temperature \t water \t CO2 \t dissolved h2o onliq \t dissolved co2 onliq \t exsolved h2o ongas \t exsolved h2o ongas mol \t gas_ontotal\n")

    # loop on pressures; volatiles calculated as initial dissolved + input flushing fluid, at 1 wt% steps
    # initial step at largest pressure
    # total initial volatile contents:
    ox[0] = previous_dissolved_h2o[0] + delta_input_gas * (1. - co2_in)
    ox[1] = previous_dissolved_co2[0] + delta_input_gas * co2_in
    
    for i in range(0, vector_size - 1):

        p = final_pressure - (delta_pressure * i)
        xy = sw.solwcad(p,T,ox,3)
        exsolved_h2o_ongas = xy[2]

        exsolved_h2o_ongas_molar = exsolved_h2o_ongas * M_co2 / (exsolved_h2o_ongas * (M_co2 - M_h2o) + M_h2o)

        gas_ontotal = (ox[0] + ox[1] - xy[0] - xy[1])/(1. - xy[0] - xy[1])
        if gas_ontotal > 1:
            gas_ontotal = 1

        with open(saveDir + '/' + saveFile + "{:.2f}".format(input_gas) + '.out','a') as fs:
            fs.write("%f \t %f \t %f \t %f \t % f \t %f \t %f \t %f \t %f \n" % (p, T, ox[0], ox[1], xy[0], xy[1], xy[2], exsolved_h2o_ongas_molar, gas_ontotal))

        ox[0] = previous_dissolved_h2o[i+1] + gas_ontotal * exsolved_h2o_ongas
        ox[1] = previous_dissolved_co2[i+1] + gas_ontotal * xy[3]

        previous_dissolved_h2o[i] = xy[0] * (1. - gas_ontotal)
        previous_dissolved_co2[i] = xy[1] * (1. - gas_ontotal)
