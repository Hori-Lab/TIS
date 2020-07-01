#!/usr/bin/env python

import sys

tk = 0.0
tc = 0.0
tc_com = 0.0
n = 0

for il, l in enumerate(open(sys.argv[1] + '.T')):
    #if il < 1000:
    #    continue
    tk += float(l.split()[1])
    tc += float(l.split()[2])
    tc_com += float(l.split()[3])
    n += 1

print(tk / float(n) , tc / float(n) , tc_com / float(n))

#nmp = int(open(sys.argv[1]+'.psf').readline().split()[0])
#
#t = 0.0
#n = 0
#for il, l in enumerate(open(sys.argv[1] + '.ts')):
#    #if il < 1000:
#    #    continue
#    if l[0] == '#' or len(l.split()) < 2:
#        continue
#    t += float(l.split()[4])
#    n += 1
#
#print(t / float(n) * (2.0 / (3.0 * nmp * 0.00198720359)))
