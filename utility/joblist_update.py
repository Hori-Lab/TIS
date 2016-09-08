#!/usr/bin/env python

import os
import shutil
import sys
import time
from filesizes import * 

'''
List of job status 
R: Running
F: job Finished and waiting for post processing
S: ready for Submission 
C: Completed
T: Terminated
X: Something wrong
'''

# Save joblist as joblist.previous
shutil.copy2('joblist', 'joblist.previous')
shutil.copy2('joblist', 'joblist.%i' % (int(time.time())))  # for saving permanently

count_R = 0 # Running
count_L = 0 # running Last run
count_F = 0 # Finished
count_P = 0 # Prepared for submission
count_S = 0 # Submitted
count_C = 0 # Completed
count_T = 0 # Terminate
count_X = 0 # something wrong

f_out = open('joblist','w')

with open('joblist.previous','r') as f:
    for l in f:
        if l[0] == '#':
            continue
        lsp = l.split() 
        jobID = int(lsp[0])
        stat  = lsp[1]
        runID = int(lsp[2])
        step  = lsp[3]

        if not stat in ('S', 'R', 'L'):
            if stat == 'C':
                f_out.write(l)
                count_C += 1
                continue
            if stat == 'X':
                f_out.write(l)
                count_X += 1
                continue
            if stat == 'T':
                f_out.write(l)
                count_T += 1
                continue
            if stat == 'F':
                f_out.write(l)
                count_F += 1
                continue
            if stat == 'P':
                f_out.write(l)
                count_P += 1
                continue
    
        path_ts  = 'f%03i/azo.ts'  % jobID
        path_dcd = 'f%03i/azo.dcd' % jobID
        path_rst = 'f%03i/azo.rst' % jobID
    
        if stat == 'S':
            if (os.path.exists(path_ts)  and 
                os.path.exists(path_dcd) and
                os.path.exists(path_rst) ):
                stat = 'R'
            else:
                count_S += 1
                f_out.write('%03i S %3i  %4s ----\n' % (jobID, runID, step))
                continue

        #if step == '----':
        #    pass
        if step == '200m':
            size_OK = size_200m
        elif runID == 0 and step == '100m':
            size_OK = size_100m_1st
        elif step == '100m':
            size_OK = size_100m
        else:
            count_X += 1
            f_out.write('%03i X %3i  %4s %3i' % (jobID, runID, step, 0) )
            f_out.write('% size is not defined\n')
            continue
    
        size_ts  = os.path.getsize(path_ts)
        size_dcd = os.path.getsize(path_dcd)
        size_rst = os.path.getsize(path_rst)

        if (size_ts  == size_OK['TS']  and 
            size_dcd == size_OK['DCD'] and
            size_rst == size_OK['RST']  ):

            # check consistency to previous run
            #with open(path_ts,'r') as fene:
            #    Rg_last = float( fene.readline().split()[-1] )
            #with open('%03i/run%i/Energies.dat' % (jobID, runID),'r') as fene:
            #    for l in fene:
            #        pass 
            #    Rg_begin = float( l.split()[-1] )
            #if not (Rg_last - 0.0001 < Rg_begin and Rg_begin < Rg_last + 0.0001):
            #    stat = 'X'
            #    count_X += 1
            #    f_out.write('%03i X %3i  %4s %3i' % (jobID, runID, step, 100*size_ts/float(size_OK['TS']) ))
            #    f_out.write('% discontinuity\n')
            #    continue

            if stat == 'L':
                count_T += 1
                f_out.write('%03i T %3i  %4s ----\n' % (jobID, runID, step) )
                continue
            elif stat == 'R':
                count_F += 1
                f_out.write('%03i F %3i  %4s ----\n' % (jobID, runID, step) )
                continue
            else:
                print ('Error: neither L nor R (1)')
                sys.exit(2)

        elif (size_ts  > size_OK['TS']  or
              size_dcd > size_OK['DCD'] or
              size_rst > size_OK['RST'] ):
            count_X += 1
            f_out.write('%03i X %3i  %4s %3i' % (jobID, runID, step, 100*size_ts/float(size_OK['TS']) ))
            f_out.write('% over-size\n')
            continue

        else:
            if stat == 'L':
                count_L += 1
                f_out.write('%03i L %3i  %4s %3i' % (jobID, runID, step, 100*size_ts/float(size_OK['TS']) ))
                f_out.write('%\n')
                continue
            elif stat == 'R':
                count_R += 1
                f_out.write('%03i R %3i  %4s %3i' % (jobID, runID, step, 100*size_ts/float(size_OK['TS']) ))
                f_out.write('%\n')
                continue
            else:
                print ('Error: neither L nor R (2)')
                sys.exit(2)

f_out.write('# Running: %i\n' % count_R)
f_out.write('# running Last run: %i\n' % count_L)
f_out.write('# Finished: %i\n' % count_F)
f_out.write('# Prepared: %i\n' % count_P)
f_out.write('# Submitted: %i\n' % count_S)
f_out.write('# Terminated: %i\n' % count_T)
f_out.write('# Completed: %i\n' % count_C)
f_out.write('# something wrong(X): %i\n' % count_X)
f_out.write('# TOTAL: %i\n' % (count_R+count_L+count_F+count_P+count_S+count_T+count_C+count_X,))
f_out.close()

print('# Running: %i' % count_R)
print('# running Last run: %i' % count_L)
print('# Finished: %i' % count_F)
print('# Prepared: %i' % count_P)
print('# Submitted: %i' % count_S)
print('# Terminated: %i' % count_T)
print('# Completed: %i' % count_C)
print('# something wrong(X): %i' % count_X)
print('# TOTAL: %i' % (count_R+count_L+count_F+count_P+count_S+count_T+count_C+count_X,))
