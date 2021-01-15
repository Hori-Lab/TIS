# Generate job directories and input files
./prepare.py


# Execute simulations
for T in 60 62 64 66 68 70 72 74 76
do
   ../../../md dimer_0.${T}/dimer.inp 1> dimer_0.${T}/log 2> dimer_0.${T}/err
done


# Analysis
rm -f dG_test.out
touch dG_test.out

for T in 60 62 64 66 68 70 72 74 76
do
   cd dimer_0.$T
   ../../st_hist.py
   head -n 1 st_hist.out >> ../dG_test.out
   cd ../
done


# Plot
gnuplot ../Figure3.gnu
