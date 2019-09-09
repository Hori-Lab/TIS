#!/usr/bin/env python

Tmin = 0.0
Tmax = 130.0

Nrep = 32

wid = (Tmax - Tmin) / float(Nrep-1)

for irep in range(1, Nrep+1):

    Tc = Tmin + (irep-1)*wid
    T = Tc + 273.15
    print('REPLICA(%i) = %5.2f' % (irep, T))
