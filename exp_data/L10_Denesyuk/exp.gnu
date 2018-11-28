set term postscript color 
set out 'exp.ps'

set size ratio 1
set key left bottom
set xr [0:140]
set yr [-0.4:-0.1]

p 'ExpChemicalPotential.txt' i 0 u 1:(0.0019872041*($1+273.15)*log($2)) title '1'\
, 'ExpChemicalPotential.txt' i 1 u 1:(0.0019872041*($1+273.15)*log($2)) title '2'\
, 'ExpChemicalPotential.txt' i 2 u 1:(0.0019872041*($1+273.15)*log($2)) title '3'\
, 'ExpChemicalPotential.txt' i 3 u 1:(0.0019872041*($1+273.15)*log($2)) title '4'\
, 'ExpChemicalPotential.txt' i 4 u 1:(0.0019872041*($1+273.15)*log($2)) title '5'\
, 'ExpChemicalPotential.txt' i 5 u 1:(0.0019872041*($1+273.15)*log($2)) title '6'
