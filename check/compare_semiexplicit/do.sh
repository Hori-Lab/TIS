cp ~/semiexplicit/examples/md.dcd semiexplicit.dcd
cp ~/semiexplicit/examples/md.out semiexplicit.out

../../md bwyv19_dcd.inp 1> log 2> err

gnuplot compare.gnu
ps2pdf compare.ps
