set term postscript color 20 lw 4 solid
set out 'Figure3.ps'

set size ratio 1

set title 'r(AG), h = 5.98, s = 5.30'
set xl 'kT / kcal/mol'
set yl 'dG / kcal/mol'

set xr [0.59:0.77]
set xtics 0.04
set mxtics 2
set yr [-1.0:1.0]
set ytics 0.2

set grid xtics ytics mxtics

unset key
p 'dG.out' u 1:3 w p lt 7 pt 7 ps 2

unset grid
set xl 'U_{ST} / kcal/mol'
set yl 'P(U_{ST})'

set xr [-7:0]
set xtics 1.0
set yr [0.0:*]
set ytics auto

set title 'kT = 0.60 kcal/mol'
p 'dimer_0.60/st_hist.out' i 1 w l lw 3

set xr [-7:0]
set title 'kT = 0.62 kcal/mol'
p 'dimer_0.62/st_hist.out' i 1 w l lw 3

set xr [-7:0]
set title 'kT = 0.64 kcal/mol'
p 'dimer_0.64/st_hist.out' i 1 w l lw 3

set xr [-7:0]
set title 'kT = 0.66 kcal/mol'
p 'dimer_0.66/st_hist.out' i 1 w l lw 3

set xr [-7:0]
set title 'kT = 0.68 kcal/mol'
p 'dimer_0.68/st_hist.out' i 1 w l lw 3

set xr [-7:0]
set title 'kT = 0.70 kcal/mol'
p 'dimer_0.70/st_hist.out' i 1 w l lw 3

set xr [-7:0]
set title 'kT = 0.72 kcal/mol'
p 'dimer_0.72/st_hist.out' i 1 w l lw 3

set xr [-7:0]
set title 'kT = 0.74 kcal/mol'
p 'dimer_0.74/st_hist.out' i 1 w l lw 3

set xr [-7:0]
set title 'kT = 0.76 kcal/mol'
p 'dimer_0.76/st_hist.out' i 1 w l lw 3
