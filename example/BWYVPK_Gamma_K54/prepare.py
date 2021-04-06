#!/usr/bin/env python
import os

conc_K = 0.054
concs_Mg = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0,]

for conc_Mg in concs_Mg:

    dirpath = 'K54_Mg%5.3f' % conc_Mg

    os.mkdir(dirpath)

    buf = 10.0   # Default
    
    if conc_Mg == 0.002:
        box = 5500.0
        cut = 200.0
        buf = 50.0

    elif conc_Mg == 0.005:
        box = 4500.0
        cut = 180.0
        buf = 40.0

    elif conc_Mg == 0.01:
        box = 3200.0
        cut = 150.0
        buf = 30.0

    elif conc_Mg == 0.02:
        box = 2600.0
        cut = 150.0
        buf = 20.0
        
    elif conc_Mg == 0.05:
        box = 1900.0
        cut = 120.0
        buf = 20.0

    elif conc_Mg == 0.1:
        box = 1500.0
        cut = 100.0

    elif conc_Mg == 0.2:
        box = 1200.0
        cut = 70.0

    elif conc_Mg == 0.5:
        box = 900.0
        cut = 70.0

    elif conc_Mg == 1.0:
        box = 700.0
        cut = 50.0

    num_Mg = int(conc_Mg * 0.001 * 6.022e-04 * box**3)
    #print (num_Mg)

    fout = open('%s/BWYVPK_Gamma.inp' % dirpath, 'w')
    for l in open('BWYVPK_Gamma_K54_temp.inp'):
        l = l.replace('##BOX##', "%f" % box)
        l = l.replace('##NUM_MG##', "%d" % num_Mg)
        l = l.replace('##BUFFER##', "%f" % buf)
        l = l.replace('##CUT_ELE##', "%f" % cut)
        fout.write(l)

    fout.close()

# \$bin/bin/runmd -p $dir/1l2x.pdb -b $dir/hbond.dat -k $dir/stack.dat -m \$bin/src/maxi_explicit -u \$bin/uvv/pmf_Mg_P_1264 -T 25 -M $mg -K $Kconc -s 5000000000 -v "$box $box $box" -a 20000 -e 20000 -C $cut -S \$RANDOM -z 29 -d $buffer $rst
