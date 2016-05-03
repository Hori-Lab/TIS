#!/usr/bin/env python

#f_out = open("azo_00_00.cg.pdb", "w")
#
#for l in open("azo_00_00.cg.pdb.original","r"):
#    if l[0:4] == 'ATOM':
#        x = float(l[30:38])
#        y = float(l[38:46])
#        z = float(l[46:54])
#        f_out.write("%s%8.3f%8.3f%8.3f%s" % (l[0:30],x-175.0,y-175.0,z-175.0,l[54:]))
#    else:
#        f_out.write(l)
#f_out.close()
#
#f_out = open("azo_05_12.ion.pdb", "w")
#for l in open("azo_05_12.ion.pdb.original","r"):
#    if l[0:4] == 'ATOM' or l[0:4] == 'HETA':
#        x = float(l[30:38])
#        y = float(l[38:46])
#        z = float(l[46:54])
#        f_out.write("%s%8.3f%8.3f%8.3f%s" % (l[0:30],x-175.0,y-175.0,z-175.0,l[54:]))
#    else:
#        f_out.write(l)
#f_out.close()

f_out = open("azo.rna.pdb", "w")
for l in open("azo_05_12.rna.pdb.original","r"):
    if l[0:4] == 'ATOM' or l[0:4] == 'HETA':
        x = float(l[30:38])
        y = float(l[38:46])
        z = float(l[46:54])
        f_out.write("%s%8.3f%8.3f%8.3f%s" % (l[0:30],x-175.0,y-175.0,z-175.0,l[54:]))
    else:
        f_out.write(l)
f_out.close()
