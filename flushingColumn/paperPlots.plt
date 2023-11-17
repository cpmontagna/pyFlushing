set ylabel 'Pressure, MPa'
set yrange [400:0] reverse

cd '../pyResults'

set term wxt 1
unset xrange
set xlabel 'Gas composition, CO_2 mole fraction'
set key bottom left
plot 'initialColumnPy.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'initial condition', 'pythonColumn0.10.out' u (1 - $8):($1*1.e-6) title 'flushing gas 10 wt%' w lines lw 2, 'pythonColumn0.50.out' u (1 - $8):($1*1.e-6) title 'flushing gas 50 wt%' w lines lw 2, 'pythonColumn1.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 100 wt%' w lines lw 2, 'pythonColumn3.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 300 wt%' w lines lw 2, 'pythonColumn5.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 500 wt%' w lines lw 2, 'pythonColumn10.00.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 1000 wt%'

set term wxt 2
set xlabel 'Mass ratio CO2/H2O in the excess fluid phase'
unset xrange
set key top right
plot 'initialColumnPy.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'initial condition', \
'pythonColumn0.10.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 10 wt%', \
'pythonColumn0.50.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 50 wt%', \
'pythonColumn1.00.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 100 wt%', \
'pythonColumn3.00.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 300 wt%', \
'pythonColumn5.00.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 500 wt%', \
     'pythonColumn10.00.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 1000 wt%'

set term wxt 3
set xlabel 'Dissolved water, wt%'
unset xrange
set key bottom left
plot 'initialColumnPy.out' using ($5*1.e2):($1*1.e-6) smooth bezier w lines lw 2 title 'initial condition', \
'pythonColumn0.10.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 10 wt%', \
'pythonColumn0.50.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 50 wt%', \
'pythonColumn1.00.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 100 wt%', \
'pythonColumn3.00.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 300 wt%', \
'pythonColumn5.00.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 500 wt%', \
     'pythonColumn10.00.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 1000 wt%'

cd '../solwcad'

set term wxt 4
unset xrange
set xlabel 'Viscosity, Pa s'
set key bottom right
plot 'initialFullColumn.out' u (10**$14):1 w lines lw 2 title 'initial condition', 'fullColumn0.10.out' u (10**$14):1 title 'flushing gas 10 wt%' w lines lw 2, 'fullColumn0.50.out' u (10**$14):1 title 'flushing gas 50 wt%' w lines lw 2, 'fullColumn1.00.out' u (10**$14):1 title 'flushing gas 100 wt%' w lines lw 2, 'fullColumn3.00.out' u (10**$14):1 title 'flushing gas 300 wt%' w lines lw 2, 'fullColumn5.00.out' u (10**$14):1 title 'flushing gas 500 wt%' w lines lw 2, 'fullColumn10.00.out' u (10**$14):1 w lines lw 2 title 'flushing gas 1000 wt%'

# comparisons

cd '../pyResults'

set term wxt 5
unset xrange
set xlabel 'Gas composition, CO_2 mole fraction'
set key bottom left
plot 'initialColumnPy.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'initial condition', 'pythonColumn0.10.out' u (1 - $8):($1*1.e-6) title 'flushing gas 10 wt%' w lines lw 2, 'pythonColumn0.50.out' u (1 - $8):($1*1.e-6) title 'flushing gas 50 wt%' w lines lw 2, 'pythonColumn1.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 100 wt%' w lines lw 2, 'pythonColumn3.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 300 wt%' w lines lw 2, 'pythonColumn5.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 500 wt%' w lines lw 2, 'pythonColumn10.00.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 1000 wt%', '../../stromboliLowCO2/pyResults/initialColumnPy.out' u (1 - $8):($1*1.e-6) title 'initial condition', '../../stromboliLowCO2/pyResults/pythonColumn0.10.out' u (1 - $8):($1*1.e-6) title 'flushing gas 10 wt%', '../../stromboliLowCO2/pyResults/pythonColumn0.50.out' u (1 - $8):($1*1.e-6) title 'flushing gas 50 wt%', '../../stromboliLowCO2/pyResults/pythonColumn1.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 100 wt%', '../../stromboliLowCO2/pyResults/pythonColumn3.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 300 wt%', '../../stromboliLowCO2/pyResults/pythonColumn5.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 500 wt%', '../../stromboliLowCO2/pyResults/pythonColumn10.00.out' u (1 - $8):($1*1.e-6) title 'flushing gas 1000 wt%'

set term wxt 6
set xlabel 'Mass ratio CO2/H2O in the excess fluid phase'
unset xrange
set key top right
plot 'initialColumnPy.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'initial condition', \
'pythonColumn0.10.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 10 wt%', \
'pythonColumn0.50.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 50 wt%', \
'pythonColumn1.00.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 100 wt%', \
'pythonColumn3.00.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 300 wt%', \
'pythonColumn5.00.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 500 wt%', \
     'pythonColumn10.00.out' using ((1. - $7)/($7)):($1*1.e-6) w lines lw 2 title 'flushing gas 1000 wt%', '../../stromboliLowCO2/pyResults/initialColumnPy.out' using ((1. - $7)/($7)):($1*1.e-6) title 'initial condition', \
'../../stromboliLowCO2/pyResults/pythonColumn0.10.out' using ((1. - $7)/($7)):($1*1.e-6) title 'flushing gas 10 wt%', \
'../../stromboliLowCO2/pyResults/pythonColumn0.50.out' using ((1. - $7)/($7)):($1*1.e-6) title 'flushing gas 50 wt%', \
'../../stromboliLowCO2/pyResults/pythonColumn1.00.out' using ((1. - $7)/($7)):($1*1.e-6) title 'flushing gas 100 wt%', \
'../../stromboliLowCO2/pyResults/pythonColumn3.00.out' using ((1. - $7)/($7)):($1*1.e-6) title 'flushing gas 300 wt%', \
'../../stromboliLowCO2/pyResults/pythonColumn5.00.out' using ((1. - $7)/($7)):($1*1.e-6) title 'flushing gas 500 wt%', \
     '../../stromboliLowCO2/pyResults/pythonColumn10.00.out' using ((1. - $7)/($7)):($1*1.e-6) title 'flushing gas 1000 wt%'

set term wxt 7
set xlabel 'Dissolved water, wt%'
unset xrange
set key bottom left
plot 'initialColumnPy.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'initial condition', \
'pythonColumn0.10.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 10 wt%', \
'pythonColumn0.50.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 50 wt%', \
'pythonColumn1.00.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 100 wt%', \
'pythonColumn3.00.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 300 wt%', \
'pythonColumn5.00.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 500 wt%', \
     'pythonColumn10.00.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 1000 wt%', '../../stromboliLowCO2/pyResults/initialColumnPy.out' using ($5*1.e2):($1*1.e-6) title 'initial condition', \
'../../stromboliLowCO2/pyResults/pythonColumn0.10.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 10 wt%', \
'../../stromboliLowCO2/pyResults/pythonColumn0.50.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 50 wt%', \
'../../stromboliLowCO2/pyResults/pythonColumn1.00.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 100 wt%', \
'../../stromboliLowCO2/pyResults/pythonColumn3.00.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 300 wt%', \
'../../stromboliLowCO2/pyResults/pythonColumn5.00.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 500 wt%', \
     '../../stromboliLowCO2/pyResults/pythonColumn10.00.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 1000 wt%'