##(1)Step     (2)E_BOND       (3)E_ANGLE      (4)E_EXCLUDED   (5)E_NON-NAT_HB (6)E_NAT_HB     (7)E_NON-NAT_ST (8)E_NAT_STACK  (9)E_ele        (10)E_potential (11)E_kinetic   (12)Rg              12:  -14
#            0   3.3239672e+01   6.7576911e+01   1.7298308e+04   0.0000000e+00  -1.1347357e+02  -3.0515079e+00  -6.3223902e+01   3.0829871e+00   1.7222458e+04   9.2979464e+01   1.2551702e+01     -5.0000
#
#u   step    tempk     radg       etot      velet  qscore     rmsd     local      go      repul      stack     tstack      hbond     thbond      elect
#      1        2       3          4          5        6      7         8          9        10         11         12        13         14         15         
#

set term postscript color 20 lw 4 solid background rgb 'white'
set out 'compare.ps'


# Total energy
set title 'Total energy'

set xr [5:*]
p 'semiexplicit.out' u 1:10 w lp title 'Hung'\
, 'md.ts' u ($1+2):4 w l title 'mine'

set xr [4000:4500]
p 'semiexplicit.out' u 1:10 w lp title 'Hung'\
, 'md.ts' u ($1+2):4 w l title 'mine'


# bond & angle
set title 'Bond and angles'

set xr [5:*]
p 'semiexplicit.out' u 1:($2+$3) w lp title 'Hung'\
, 'md.ts' u ($1+2):8 w l title 'mine'

set xr [4000:4500]
p 'semiexplicit.out' u 1:($2+$3) w lp title 'Hung'\
, 'md.ts' u ($1+2):8 w l title 'mine'


# H-bonds
set title 'H-bonds'

set xr [5:*]
p 'semiexplicit.out' u 1:($5+$6) w lp title 'Hung'\
, 'md.ts' u ($1+2):($13+$14) w l title 'mine'

set xr [4000:4500]
p 'semiexplicit.out' u 1:($5+$6) w lp title 'Hung'\
, 'md.ts' u ($1+2):($13+$14) w l title 'mine'


# Stack
set title 'Stack'

set xr [5:*]
p 'semiexplicit.out' u 1:($7+$8) w lp title 'Hung'\
, 'md.ts' u ($1+2):($11+$12) w l title 'mine'

set xr [4000:4500]
p 'semiexplicit.out' u 1:($7+$8) w lp title 'Hung'\
, 'md.ts' u ($1+2):($11+$12) w l title 'mine'


# Excluded volume
set title 'Excluded volume'

set xr [5:*]
p 'semiexplicit.out' u 1:4 w lp title 'Hung'\
, 'md.ts' u ($1+2):10 w l title 'mine'

set xr [4000:4500]
p 'semiexplicit.out' u 1:4 w lp title 'Hung'\
, 'md.ts' u ($1+2):10 w l title 'mine'


# Electrostaic
set title 'Electrostatic'

set xr [5:*]
p 'semiexplicit.out' u 1:9 w lp title 'Hung'\
, 'md.ts' u ($1+2):15 w l title 'mine'

set xr [4000:4500]
p 'semiexplicit.out' u 1:9 w lp title 'Hung'\
, 'md.ts' u ($1+2):15 w l title 'mine'
