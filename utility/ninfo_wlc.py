#!/usr/bin/env python

N = 100

for i in range(1,N):
    print ('fene %5i  1  1 %3i %3i %3i %3i    3.8   2.0   20.0000' % (i, i, i+1, i, i+1))

for i in range(1,N-1):
    print ('angl %5i  1  1 %3i %3i %3i %3i %3i %3i    170.0000  1.0  1.0     5.0000' % (i, i, i+1, i+2, i, i+1, i+2))
