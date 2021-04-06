set term postscript color 20 lw 4 solid
set out 'BWYVPK_Gamma_K54.ps'

set logscale x
set grid

set xl 'C_{bulk} / mM'
set yl 'Gamma'

set xr [0.001:1.2]

unset key
p 'gamma.out' u ($1*1000):2 w p pt 7 ps 2

