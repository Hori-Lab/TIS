set term postscript color 20 lw 4 solid
set out 'fit.ps'
set yr [-5:12]
set xr [1:10]
set xl 'r'
set yl 'W(r)'
p 'fit_2_4exp_sig.out' u 1:2 w l lw 2 title 'W(r)'\
, 'fit_2_4exp_sig.out' u 1:3 w lp title 'Fit'

set out 'fit_error.ps'
set yl 'Fit(r) - W(r)'
set yr [-0.1:0.1]
p 'fit_2_4exp_sig.out' u 1:($3-$2) w l lw 1 title 'Error'

