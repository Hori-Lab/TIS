~/Futil/TIS_analysis/ion_gamma K54_Mg0.002/md.dcd 50000 85 200 5500.0 > K54_Mg0.002/gamma.out
~/Futil/TIS_analysis/ion_gamma K54_Mg0.005/md.dcd 50000 85 274 4500.0 > K54_Mg0.005/gamma.out
~/Futil/TIS_analysis/ion_gamma K54_Mg0.010/md.dcd 50000 85 197 3200.0 > K54_Mg0.010/gamma.out
~/Futil/TIS_analysis/ion_gamma K54_Mg0.020/md.dcd 50000 85 211 2600.0 > K54_Mg0.020/gamma.out
~/Futil/TIS_analysis/ion_gamma K54_Mg0.050/md.dcd 50000 85 206 1900.0 > K54_Mg0.050/gamma.out
~/Futil/TIS_analysis/ion_gamma K54_Mg0.100/md.dcd 50000 85 203 1500.0 > K54_Mg0.100/gamma.out
~/Futil/TIS_analysis/ion_gamma K54_Mg0.200/md.dcd 50000 85 208 1200.0 > K54_Mg0.200/gamma.out
~/Futil/TIS_analysis/ion_gamma K54_Mg0.500/md.dcd 50000 85 219  900.0 > K54_Mg0.500/gamma.out
~/Futil/TIS_analysis/ion_gamma K54_Mg1.000/md.dcd 50000 85 206  700.0 > K54_Mg1.000/gamma.out

tail -qn1 K54_Mg*/gamma.out > gamma.out
