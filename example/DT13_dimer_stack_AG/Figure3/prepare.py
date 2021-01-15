#!/usr/bin/env python

import os

for kT in (0.6, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76):
    T = kT / 0.00198720359 

    dirname = './dimer_%4.2f' % kT
    os.mkdir(dirname)

    f_inp = open('%s/dimer.inp' % (dirname,), 'w')
    for l in open('dimer_template.inp'):
        l = l.replace('##DIR##', dirname).replace('##TEMP##', '%5.2f' % (T,))
        f_inp.write(l)
    
    f_inp.close()

