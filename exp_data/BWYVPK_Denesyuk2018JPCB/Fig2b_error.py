#!/usr/bin/env python

import math

def simulation_line(logC):
    return -4.91268 * logC - 16.92273

data_logC = []
data_dG = []

s2 = 0.0
n = 0
for l in open('BWYVPK_Fig2b.txt'):
    lsp = l.split()
    logC = float(lsp[0])
    dG_exp = float(lsp[1])

    dG_sim = simulation_line(logC)

    err = dG_exp - dG_sim
    print err 

    s2 += err ** 2
    n += 1

print 'MSE = ', s2 / float(n)
print 'RMSD = ',  math.sqrt(s2 / float(n))
