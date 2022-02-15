set ylabel 'Pressure, MPa'
set yrange [405:0.09] reverse

set term qt 1
unset xrange
set xlabel 'Gas composition, CO_2 mole fraction'
set key bottom left
plot 'volumes/initialColumn.out' u ($8*1e-2):1 w lines lw 2 title 'initial condition', 'pythonColumn0.1.out' u (1 - $8):($1*1.e-6) title 'flushing gas 10 wt%' w lines lw 2, 'pythonColumn0.5.out' u (1 - $8):($1*1.e-6) title 'flushing gas 50 wt%' w lines lw 2, 'pythonColumn1.0.out' u (1 - $8):($1*1.e-6) title 'flushing gas 100 wt%' w lines lw 2, 'pythonColumn3.0.out' u (1 - $8):($1*1.e-6) title 'flushing gas 300 wt%' w lines lw 2, 'pythonColumn5.0.out' u (1 - $8):($1*1.e-6) title 'flushing gas 500 wt%' w lines lw 2, 'pythonColumn10.0.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 1000 wt%'

set term qt 2
set xlabel 'Dissolved water, wt%'
unset xrange
set key bottom left
plot 'volumes/initialColumn.out' using 5:1 smooth bezier w lines lw 2 title 'initial condition', \
'pythonColumn0.1.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 10 wt%', \
'pythonColumn0.5.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 50 wt%', \
'pythonColumn1.0.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 100 wt%', \
'pythonColumn3.0.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 300 wt%', \
'pythonColumn5.0.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 500 wt%', \
'pythonColumn10.0.out' using ($5*1.e2):($1*1.e-6) w lines lw 2 title 'flushing gas 1000 wt%'