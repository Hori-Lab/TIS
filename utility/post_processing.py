#!/usr/bin/env python

'''
Process jobs finised ("F" status), and prepare for next run (to be ready for submisstion "S").
'''

import shutil
import os
import sys
import time
from filesizes import *

# Save joblist as joblist.previous
shutil.copy2('joblist','joblist.previous')
shutil.copy2('joblist', 'joblist.%i' % (int(time.time())))  # for saving permanently

count_R = 0 # Running
count_L = 0 # running Last run
count_F = 0 # Finished
count_P = 0 # Prepared
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

        if not stat in ('F','T'):
            f_out.write(l)
            if stat == 'R':
                count_R += 1 
            elif stat == 'L':
                count_L += 1
            elif stat == 'S':
                count_S += 1 
            elif stat == 'P':
                count_P += 1
            elif stat == 'C':
                count_C += 1
            elif stat == 'X':
                count_X += 1
            continue

        path_ts  = 'f%03i/azo.ts' % jobID
        path_dcd = 'f%03i/azo.dcd' % jobID
        path_rst = 'f%03i/azo.rst' % jobID

        if not (os.path.exists(path_ts)  and 
                os.path.exists(path_dcd) and
                os.path.exists(path_rst) ):
            count_X += 1
            f_out.write('%03i X %3i  %4s' % (jobID, runID, step) )
            f_out.write('% missing file\n')
            continue

        if step == '----':
            pass
        elif step == '200m':
            size_OK = size_200m
        elif runID == 0 and step == '100m':
            size_OK = size_100m_1st
        elif step == '100m':
            size_OK = size_100m
        else:
            count_X += 1
            f_out.write('%03i X %3i  %4s' % (jobID, runID, step) )
            f_out.write(' size is not defined\n')
            continue
    
        size_ts  = os.path.getsize(path_ts)
        size_dcd = os.path.getsize(path_dcd)
        size_rst = os.path.getsize(path_rst)

        if (size_ts  != size_OK['TS']  or
            size_dcd != size_OK['DCD'] or
            size_rst != size_OK['RST'] ):
            count_X += 1
            f_out.write('%03i X %3i  %4s' % (jobID, runID, step) )
            f_out.write('% wrong file size\n')
            continue

        # Generate next run directory
        newdir = 'f%03i/run%03i' % (jobID, runID+1)
        if (os.path.exists(newdir)):
            count_X += 1
            f_out.write('%03i X %3i  %4s' % (jobID, runID, step) )
            f_out.write(' next run directory exists\n')
            continue

        os.mkdir(newdir)

        # Move data to the new directory
        shutil.move(path_ts, newdir)
        shutil.move(path_dcd, newdir)
        shutil.move(path_rst, newdir)
        shutil.move('f%03i/azo.data' % (jobID,), newdir)
        shutil.move('f%03i/azo.hb'   % (jobID,), newdir)
        shutil.move('f%03i/azo.st'   % (jobID,), newdir)
        shutil.move('f%03i/azo.ninfo'% (jobID,), newdir)
        shutil.move('f%03i/azo.psf'  % (jobID,), newdir)
        shutil.move('f%03i/azo.pdb'  % (jobID,), newdir)
        shutil.move('f%03i/f%03i.err' % (jobID, jobID), newdir)
        shutil.move('f%03i/f%03i.log' % (jobID, jobID), newdir)
        shutil.move('f%03i/f%03i.inp' % (jobID, jobID), newdir)
        if (os.path.exists('f%03i/azo.restart' % (jobID,))):
            shutil.move('f%03i/azo.restart' % (jobID,), newdir)   # This is just a link
        if (os.path.exists('f%03i/%03i.rst' % (jobID,runID))):
            shutil.move('f%03i/%03i.rst'  % (jobID, runID), newdir)

        if stat == 'T':
            count_C += 1
            f_out.write('%03i C %3i  ---- ----\n' % (jobID, runID+1, ) )
            continue

        elif stat == 'F':
            # Copy rst file
            shutil.copy2('f%03i/run%03i/azo.rst' % (jobID, runID+1), 'f%03i/%03i.rst' % (jobID, runID+1))
            os.symlink('%03i.rst' % (runID+1,), 'f%03i/azo.restart' % (jobID,))

            # Generate new inp
            f_newinp = open('f%03i/f%03i.inp' % (jobID, jobID), 'w')
            with open('f%03i/run%03i/f%03i.inp' % (jobID, runID+1, jobID), 'r') as f_preinp:
                for l_pre in f_preinp:
                    if l_pre.find('n_tstep(1)') != -1:
                        f_newinp.write('n_tstep(1) = %i\n' % ((runID+1)*200000000))
                    else:
                        f_newinp.write(l_pre)

            count_P += 1
            if runID+1 < 2:
                f_out.write('%03i P %3i  %4s ----\n' % (jobID, runID+1, '100m') )
            else:
                f_out.write('%03i P %3i  %4s ----\n' % (jobID, runID+1, '200m') )
            continue

        else:
            print ('Error: neither T nor F')
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
