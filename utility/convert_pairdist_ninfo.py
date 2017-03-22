#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print 'Usage: SCRIPT [pairdist file]'
    sys.exit(2)

i_id = 0
for l in open(sys.argv[1],'r'):
    lsp = l.split()
    i = int(lsp[0])
    j = int(lsp[1])
    d = float(lsp[2])
    if abs(i-j) > 2 and d < 8.0:
        i_id += 1
        #print 'con_gauss  %i 1 1 %i %i %i %i' % (i_id, i,j,i,j)
        print 'wca  %i 1 1 %i %i %i %i  6.30  0.59645916047  0.59645916047' % (i_id, i,j,i,j)

