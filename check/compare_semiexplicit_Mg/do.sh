# Run the following in semiexplicit
#$bin/semiexplicit -p 1l2x.pdb -b hbond.dat -k stack.dat -m $bin/maxi_explicit -u $bin/uvv/pmf_Mg_P_1264 -T 25 -M 0.001 -K 0.057 -E 1000000 -s 10000 -v "700 700 700" -a 1 -e 1 -C 50 -S 5 -z 29 -d 10 -P bwyv_K57mM_Mg1mM.pdb -o bwyv_K57mM_Mg1mM.out -t bwyv_K57mM_Mg1mM.dcd -i bwyv_K57mM_Mg1mM.mdinfo -r bwyv_K57mM_Mg1mM.rst 1> bwyv_K57mM_Mg1mM.log 2> bwyv_K57mM_Mg1mM.err

cp ~/semiexplicit/examples/bwyv_K57mM_Mg1mM.dcd semiexplicit.dcd
cp ~/semiexplicit/examples/bwyv_K57mM_Mg1mM.out semiexplicit.out

../../md bwyv19_Mg_dcd.inp 1> log 2> err

gnuplot compare.gnu
ps2pdf compare.ps
