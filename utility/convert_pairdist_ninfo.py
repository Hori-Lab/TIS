#!/usr/bin/env python

i_id = 0
for l in open('1vii.a.pairdist.dat','r'):
    lsp = l.split()
    i = int(lsp[0])
    j = int(lsp[1])
    d = float(lsp[2])
    if abs(i-j) > 2 and d <= 8.5:
        i_id += 1
        print 'con_gauss  %i 1 1 %i %i %i %i' % (i_id, i,j,i,j)

