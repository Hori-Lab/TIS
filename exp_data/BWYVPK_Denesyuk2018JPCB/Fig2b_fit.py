#!/usr/bin/env python

import scipy.stats

data_logC = []
data_dG = []

for l in open('BWYVPK_Fig2b.txt'):
    lsp = l.split()
    data_logC.append( float(lsp[0]) )
    data_dG.append(   float(lsp[1]) )

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(data_logC, data_dG)

print 'slope', slope
print 'intercept', intercept
print 'r_value', r_value
print 'p_value',  p_value
print 'std_err=', std_err
