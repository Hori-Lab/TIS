#!/usr/bin/env python

N_AVO = 6.022140857e23   #< Avogadro constant [/mol]


# Fix the number of Mg
concMg = 0.002
numMg = 220
volume = float(numMg) / (concMg * (N_AVO * 1.0e-27))
box = volume ** (1./3.)
print (box)


# Fix the box size
concMg = 0.002
box = 550.0
volume = box ** 3
numMg = concMg * (N_AVO * 1.0e-27) * volume
print (numMg)

