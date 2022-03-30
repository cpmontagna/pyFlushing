import solwcadsub as sw        # import solwcad routine
import numpy as np
import os

# ============ INPUTS =========================================================

T = 1473.                      # temperature, K
initial_pressure = 1.e5
final_pressure = 4.e8
delta_pressure = 1.e7

# magma composition: weight fractions of water, CO2 and ten major oxides (SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O)
# silicate melt composition: sample ST531 4b Metrich et al. 2010
# initial water and CO2 have no role as they are computed from equilibrium with a specific fluid phase composition
oxidesList = [.02, .005, .4741, .0087, .1496, .0164, .0532, .0009, .0394, .1556, .022, .016]

# composition (CO2 wt frac) of the fluid phase in equilibrium with the initial column
wCO2 = 0.85
# file to read for equilibrium data
inFile = 'solwcadResultsIsopleth.out'

# flushing phase
# initial wt fraction of input flushing fluid
initial_input_gas = 0.01
final_input_gas = 10.
delta_input_gas = 0.01 
co2_in = 1.0

# save results
saveDir = 'equilibrium4'
initialCondition = 'equilibriumColumnPy.out' 
saveFile = 'equilibriumColumn'

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
equilibriumIndex = np.zeros(vector_size)

# initial condition: magma column in equilibrium with a fluid phase of given composition
with open(saveDir + '/' + initialCondition,'w') as fs:
    fs.write("weight fractions\n")
    fs.write("pressure \t dissolvedH2Oonliq \t dissolvedCO2onliq \t exsolvedCO2ongas \n")

# load data from file
eq_p,eq_diss_h2o,eq_diss_co2,eq_exs_co2 = np.loadtxt(inFile, usecols = (0,4,5,7), unpack = True)
eq_p_unique, indeces = np.unique(eq_p, return_index = True)

# check consistency in size of pressure vectors
if eq_p_unique.size != vector_size:
    print('check pressure vectors')

indeces = np.sort(indeces)
pStep = indeces[1] - indeces[0]

# find nearest wCO2 for each pressure
i = 0
for j in range(0, eq_p.size, pStep):
    wCO2step = eq_exs_co2[j:j + pStep]
    [ind, minCO2] = min(enumerate(wCO2step), key=lambda x: abs(x[1] - wCO2))
    equilibriumIndex[i] = j + ind
    i = i+1

for i in range(0,vector_size):
    p = eq_p[int(equilibriumIndex[i])]
    previous_dissolved_h2o[i] = eq_diss_h2o[int(equilibriumIndex[i])]
    previous_dissolved_co2[i] = eq_diss_co2[int(equilibriumIndex[i])]
    initial_exs_co2 = eq_exs_co2[int(equilibriumIndex[i])]    
        
    with open(saveDir + '/' + initialCondition,'a') as fs:
        fs.write("%f \t %f \t %f \t %f \n" % (p, previous_dissolved_h2o[i], previous_dissolved_co2[i], initial_exs_co2))

# flushing: at each pressure step, the excess fluid phase is added to
# dissolved volatiles in the next step
for j in range(0, input_size):
    
    input_gas = initial_input_gas + delta_input_gas * j
    with open(saveDir + '/' + saveFile + str(input_gas) + '.out','w') as fs:
        fs.write("weight fractions unless otherwise noted\n")
        fs.write("pressure \t temperature \t water \t CO2 \t dissolvedH2Oonliq \t dissolvedCO2onliq \t exsolvedH2Oongas \t exsolvedH2OongasMol \t gasOntot\n")

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

        with open(saveDir + '/' + saveFile + str(input_gas) + '.out','a') as fs:
            fs.write("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n" % (p, T, ox[0], ox[1], xy[0], xy[1], xy[2], exsolved_h2o_ongas_molar, gas_ontotal))

        ox[0] = previous_dissolved_h2o[i+1] + gas_ontotal * exsolved_h2o_ongas
        ox[1] = previous_dissolved_co2[i+1] + gas_ontotal * xy[3]

        previous_dissolved_h2o[i] = xy[0] * (1. - gas_ontotal)
        previous_dissolved_co2[i] = xy[1] * (1. - gas_ontotal)

