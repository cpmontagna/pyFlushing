cd '../pyResults'

set ylabel 'pressure (MPa)'
set yrange [405:0.09] reverse

set term wxt 1
unset xrange
set xlabel 'total water, wt%'
set key top right
plot 'initialColumnPy.out' u ($3*1.e2):($1*1.e-6) title 'initial column'

set term wxt 2
unset xrange
set xlabel 'total water, wt%'
set key top right
plot 'pythonColumn0.01.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 1 wt%', \
'pythonColumn0.02.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 2 wt%', \
'pythonColumn0.03.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 3 wt%', \
'pythonColumn0.05.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 5 wt%', \
'pythonColumn0.11.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 10 wt%', \
'pythonColumn0.2.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 20 wt%', \
'pythonColumn0.5.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 50 wt%', \
'pythonColumn0.9.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 90 wt%', \
'pythonColumn1.1.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 110 wt%', \
'pythonColumn1.5.out' using ($3*1.e2):($1*1.e-6) title 'flushing gas 150 wt%'

set term wxt 3
unset xrange
set xlabel 'total CO2, wt%'
set key top right
plot 'initialColumnPy.out' u ($4*1.e2):($1*1.e-6) title 'initial column'

set term wxt 4
unset xrange
set xlabel 'total CO2, wt%'
set key top right
plot 'pythonColumn0.01.out' using ($4*1.e2):($1*1.e-6) title 'flushing gas 1 wt%', \
'pythonColumn0.02.out' using ($4*1.e2):($1*1.e-6) title 'flushing gas 2 wt%', \
'pythonColumn0.03.out' using ($4*1.e2):($1*1.e-6) title 'flushing gas 3 wt%', \
'pythonColumn0.05.out' using ($4*1.e2):($1*1.e-6) title 'flushing gas 5 wt%', \
'pythonColumn0.11.out' using ($4*1.e2):($1*1.e-6) title 'flushing gas 10 wt%', \
'pythonColumn0.2.out' using ($4*1.e2):($1*1.e-6) title 'flushing gas 20 wt%', \
'pythonColumn0.5.out' using ($4*1.e2):($1*1.e-6) title 'flushing gas 50 wt%', \
'pythonColumn0.9.out' using ($4*1.e2):($1*1.e-6) title 'flushing gas 90 wt%'

set term wxt 5
unset xrange
set xlabel 'dissolved water, wt%'
set key top right
plot 'initialColumnPy.out' u ($5*1.e2):($1*1.e-6) title 'initial column'

set term wxt 6
unset xrange
set xlabel 'dissolved water, wt%'
set key bottom left
plot 'pythonColumn0.01.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 1 wt%', \
'pythonColumn0.02.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 2 wt%', \
'pythonColumn0.11.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 10 wt%', \
'pythonColumn0.9.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 90 wt%', \
'pythonColumn1.1.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 110 wt%', \
'pythonColumn2.0.out' using ($5*1.e2):($1*1.e-6) title 'flushing gas 200 wt%'

set term wxt 7
unset xrange
set xlabel 'dissolved CO2, wt%'
set key top right
plot 'initialColumnPy.out' u ($6*1.e2):($1*1.e-6) title 'initial column'

set term wxt 8
unset xrange
set xlabel 'dissolved CO2, wt%'
set key top right
plot 'pythonColumn0.01.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 1 wt%', \
'pythonColumn0.02.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 2 wt%', \
'pythonColumn0.03.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 3 wt%', \
'pythonColumn0.05.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 5 wt%', \
'pythonColumn0.11.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 10 wt%', \
'pythonColumn0.2.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 20 wt%', \
'pythonColumn0.5.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 50 wt%', \
'pythonColumn0.9.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 90 wt%', \
'pythonColumn1.1.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 110 wt%', \
'pythonColumn1.5.out' using ($6*1.e2):($1*1.e-6) title 'flushing gas 150 wt%'

set term wxt 9
set xlabel 'gas composition, water wt%'
plot 'pythonColumn0.01.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 1 wt%', \
'pythonColumn0.02.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 2 wt%', \
'pythonColumn0.03.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 3 wt%', \
'pythonColumn0.05.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 5 wt%', \
'pythonColumn0.11.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 10 wt%', \
'pythonColumn0.2.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 20 wt%', \
'pythonColumn0.5.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 50 wt%', \
'pythonColumn0.9.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 90 wt%', \
'pythonColumn1.1.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 110 wt%', \
'pythonColumn1.5.out' using ($7*1.e2):($1*1.e-6) title 'flushing gas 150 wt%'

set term wxt 10
unset xrange
set xlabel 'gas composition, water mole fraction'
set key bottom right
plot 'pythonColumn0.01.out' u 8:($1*1.e-6) title 'flushing gas 1 wt%', \
'pythonColumn0.02.out' u 8:($1*1.e-6) title 'flushing gas 2 wt%', \
'pythonColumn0.03.out' u 8:($1*1.e-6) title 'flushing gas 3 wt%', \
'pythonColumn0.05.out' u 8:($1*1.e-6) title 'flushing gas 5 wt%', \
'pythonColumn0.11.out' u 8:($1*1.e-6) title 'flushing gas 10 wt%', \
'pythonColumn0.2.out' u 8:($1*1.e-6) title 'flushing gas 20 wt%', \
'pythonColumn0.5.out' u 8:($1*1.e-6) title 'flushing gas 50 wt%', \
'pythonColumn0.9.out' u 8:($1*1.e-6) title 'flushing gas 90 wt%', \
'pythonColumn1.0.out' u 8:($1*1.e-6) title 'flushing gas 100 wt%', \
'pythonColumn1.1.out' u 8:($1*1.e-6) title 'flushing gas 110 wt%', \
'pythonColumn1.5.out' u 8:($1*1.e-6) title 'flushing gas 150 wt%'

set term wxt 11
unset xrange
set xlabel 'gas composition, CO2 mole fraction'
plot 'pythonColumn0.01.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 1 wt%', \
'pythonColumn0.02.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 2 wt%', \
'pythonColumn0.03.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 3 wt%', \
'pythonColumn0.05.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 5 wt%', \
'pythonColumn0.1.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 10 wt%', \
'pythonColumn0.5.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 20 wt%', \
'pythonColumn0.5.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 50 wt%', \
'pythonColumn0.9.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 90 wt%', \
'pythonColumn1.1.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 110 wt%', \
'pythonColumn1.5.out' u (1 - $8):($1*1.e-6) w lines lw 2 title 'flushing gas 150 wt%'

set term qt 0
unset xrange
set xlabel 'Gas on total, wt%'
plot 'initialColumnPy.out' using ($9*1.e2):($1*1.e-6) w lines lw 2 title '', 'pythonColumn0.1.out' using ($9*1.e2):($1*1.e-6) w lines lw 2 title '', 'pythonColumn0.5.out' using ($9*1.e2):($1*1.e-6) w lines lw 2 title '', 'pythonColumn1.0.out' using ($9*1.e2):($1*1.e-6) w lines lw 2 title '', 'pythonColumn3.0.out' using ($9*1.e2):($1*1.e-6) w lines lw 2 title '', 'pythonColumn5.0.out' using ($9*1.e2):($1*1.e-6) w lines lw 2 title '', 'pythonColumn10.0.out' using ($9*1.e2):($1*1.e-6) w lines lw 2 title ''

set term wxt 13
unset xrange
set xlabel 'fluid CO2, wt%'
set key top right
plot 'pythonColumn0.01.out' using (($4 - $6)*1.e2):($1*1.e-6) title 'flushing gas 1 wt%', \
'pythonColumn0.02.out' using (($4 - $6)*1.e2):($1*1.e-6) title 'flushing gas 2 wt%', \
'pythonColumn0.03.out' using (($4 - $6)*1.e2):($1*1.e-6) title 'flushing gas 3 wt%', \
'pythonColumn0.05.out' using (($4 - $6)*1.e2):($1*1.e-6) title 'flushing gas 5 wt%', \
'pythonColumn0.11.out' using (($4 - $6)*1.e2):($1*1.e-6) title 'flushing gas 10 wt%', \
'pythonColumn0.2.out' using (($4 - $6)*1.e2):($1*1.e-6) title 'flushing gas 20 wt%', \
'pythonColumn0.5.out' using (($4 - $6)*1.e2):($1*1.e-6) title 'flushing gas 50 wt%', \
'pythonColumn0.9.out' using (($4 - $6)*1.e2):($1*1.e-6) title 'flushing gas 90 wt%'

set term qt 0
set xlabel 'Dissolved water, wt%'
set xrange [0:8]
set ylabel 'Dissolved CO_2, ppm'
set yrange [0:7000]
plot 'initialColumnPy.out' u ($5*1.e2):($6*1.e6)  w lines lw 2 title 'initial condition', 'pythonColumn0.1.out' u ($5*1.e2):($6*1.e6)  w lines lw 2 title 'flushing gas 10 wt%', 'pythonColumn0.5.out' u ($5*1.e2):($6*1.e6)  w lines lw 2 title 'flushing gas 50 wt%', 'pythonColumn1.0.out' u ($5*1.e2):($6*1.e6)  w lines lw 2 title 'flushing gas 100 wt%', 'pythonColumn3.0.out' u ($5*1.e2):($6*1.e6)  w lines lw 2 title 'flushing gas 300 wt%', 'pythonColumn5.0.out' u ($5*1.e2):($6*1.e6)  w lines lw 2 title 'flushing gas 500 wt%', 'pythonColumn10.0.out' u ($5*1.e2):($6*1.e6)  w lines lw 2 title 'flushing gas 1000 wt%'