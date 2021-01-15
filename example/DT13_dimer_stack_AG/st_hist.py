#!/usr/bin/env python

import numpy as np
import math

frame_skip = 1000
f_out = open('st_hist.out', 'w')


### Read T form inp file
T = 0.0

for l in open('dimer.inp'):
    if l.startswith('tempk'):
        T = float(l.split()[2])

kT = 0.00198720359 * T


### Read stall file
data = []

iframe = 0
for l in open('md.stall'):
    iframe += 1
    if iframe < frame_skip:
        continue
    
    data.append(float(l))


# dG
n = 0
for e in data:
    if e < -kT:
        n += 1
p = n / len(data)

dG = -kT * math.log(p) + kT * math.log(1.0 - p)

f_out.write('%4.2f %6.2f %f\n' % (kT, T, dG))
f_out.write('\n\n')




# Histogram
# -10 to 0
bins = [-10 + 0.1*i for i in range(0,100)]
bins.append(0.0)


H, B = np.histogram(data, bins, density=True)

for i in range(len(H)):
    f_out.write('%f %f\n' % (0.5*(B[i]+B[i+1]), H[i]))

f_out.close()
