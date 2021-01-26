#!/usr/bin/env python

ss_D = 0.90

#BOLTZ_J = 1.38064852e-23   #!< Boltzmann constant [J/K]
#N_AVO   = 6.022140857e23   #!< Avogadro constant [/mol]
#KCAL2JOUL = 4184.0         #!< (kcal -> J)  [J/kcal]
#JOUL2KCAL = 1.0/KCAL2JOUL  #!< (J -> kcal)  [kcal/J]
#JOUL2KCAL_MOL  = JOUL2KCAL * N_AVO  #!< (J -> kcal/mol)
#BOLTZ_KCAL_MOL = BOLTZ_J * JOUL2KCAL_MOL   #!< Boltzmann constant [kcal/mol/K]

#KELVIN_TO_KT = 
BOLTZ_KCAL_MOL = 1.9872057946353295439632e-3

NN = 'AA'
s = -0.319
h = 5.194 - ss_D / 0.7093838769
Tm = 0.594 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'AC'                                
s = -0.319
h = 5.146 - ss_D / 0.7185955600
Tm = 0.594 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'AG'                                  
s = 5.301
h = 5.977 - ss_D / 0.6968019829
Tm = 0.678 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'AU'                                  
s = -0.319
h = 5.146 - ss_D / 0.7185955600
Tm = 0.594 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'CA'                                  
s = -0.319
h = 5.163 - ss_D / 0.6847830171
Tm = 0.594 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'CC'                                  
s = -1.567
h = 4.873 - ss_D / 0.6991615586
Tm = 0.568 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'CG'                                  
s = 0.774
h = 5.482 - ss_D / 0.6816268897
Tm = 0.627 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'CU'                                  
s = -1.567
h = 4.873 - ss_D / 0.6832570771
Tm = 0.568 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'GA'                                  
s = 5.301
h = 5.948 - ss_D / 0.6903176657
Tm = 0.678 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'GC'                                  
s = 4.370
h = 5.927 - ss_D / 0.7042060343
Tm = 0.682 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'GG'                                  
s = 7.346
h = 6.416 - ss_D / 0.6971421514
Tm = 0.728 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'GU'                                  
s = 2.924
h = 5.836 - ss_D / 0.6984426543
Tm = 0.672 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'UA'                                  
s = -0.319
h = 5.163 - ss_D / 0.6847830171
Tm = 0.594 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'UC'                                
s = -1.567
h = 4.880 - ss_D / 0.6758595771
Tm = 0.568 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'UG'                                  
s = 2.924
h = 5.886 - ss_D / 0.7025528229
Tm = 0.672 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))

NN = 'UU'                                  
s = -3.563
h = 4.267 - ss_D / 0.6686014771
Tm = 0.500 / BOLTZ_KCAL_MOL
print ('%2s  %7.5f  %6.3f  %7.3f' % (NN, h,s,Tm))
