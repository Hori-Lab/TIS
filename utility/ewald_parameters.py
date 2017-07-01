#!/usr/bin/env python

import subprocess32

max_time = 100.0
inp_tmpl    = 'NaCl50_ewld_energy.inp'
inp_tmpl_md = 'NaCl50_ewld_energy_md.inp'

exe1 = '../md.ewld'
exe2 = '../md'

accuracy = 0.001  # kcal/mol
accuracy_for_exact = 0.000001  # kcal/mol
max_cut = 75.0  # Half of the box size


f_log = open('ewald_param.out', 'w')


'''
Find a converged energy (Econv) with alpha fixed to "middle_alpha"
The value of hmax is increased until energy value is converged with the accuracy (accuracy_for_exact)
'''
middle_alpha = 0.3
hmax_series = [1.0,2.0,5.0,10.0, 12.5,15.0,17.5, 20.0,25.0,30.0,35.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,125.0,150.0,175.0,200.0]
f_log.write('##### Find the converged energy value\n')
f_log.write('## set accuracy threshold to %f\n' % accuracy_for_exact)
f_log.write('## set cutoff_ele to %f\n' % max_cut)
f_log.write('## set alpha to %f\n' % middle_alpha)
f_log.write('#  Hmax    Energy\n')
Epre = 0.0
for hmax in hmax_series:

    ''' Prepare an input file '''
    f_inp = open('tmp.inp', 'w')
    for  l in open(inp_tmpl):
        if l.find('cutoff_ele') != -1:
            f_inp.write('cutoff_ele = %10.5f\n' % max_cut)
        elif l.find('ewld_alpha') != -1:
            f_inp.write('ewld_alpha = %10.5f\n' % middle_alpha)
        elif l.find('ewld_hmax') != -1:
            f_inp.write('ewld_hmax = %10.5f\n' % hmax)
        else:
            f_inp.write(l)
    f_inp.close()

    ''' Execute the code '''
    f_out = open('out','w')
    f_err = open('err','w')
    cmds = [ exe1, 'tmp.inp' ]
    subprocess32.call( cmds, stdout=f_out, stderr=f_err)
    f_out.close()
    f_err.close()

    ''' Read the Ewald energy output '''
    for l in open('out'):
        if l.find("EWLD: ") != -1:
            energy = float( l.strip().split()[-1] )
            break

    f_log.write('%7.3f  %13.8f\n' % (hmax, energy))
    f_log.flush()

    ''' Test if it is converged '''
    if abs(energy - Epre) <= accuracy_for_exact:
        Econv = energy
        break
    else:
        Epre = energy

f_log.write('# Econv (converged energy) = %13.8f\n' % Econv)
f_log.write('\n\n')
f_log.flush()

#Econv = -24.9563800623


'''
Test a set of combinations of (cutoff_ele, alpha, hmax) if they are accurate enough
'''
data = []
f_log.write('##### Test each combination of (alpha, hmax) whether it is accurate enough.\n')
f_log.write('## set accuracy threshold to %f\n' % accuracy_for_exact)
f_log.write('## set cutoff_ele to %f\n' % max_cut)
hmax_series = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.5, 25.0,27.5, 30.0, 40.0, 50.0, 75.0, 100.0]
alpha_series = []
alpha_max = 0.20
i = 1 ; d_alpha = 0.02
alpha = d_alpha
while (alpha <= alpha_max):
    alpha_series.append(alpha)
    i += 1
    alpha = d_alpha * i
#print '# Alpha: ', alpha_series

f_log.write('# alpha     Hmax      Energy       abs(Energy-Econv) T(main_loop)\n')
for alpha in alpha_series:
    flg_skip = False
    Epre = 0.0
    for hmax in hmax_series:

        if flg_skip:
            data.append( (alpha, hmax, Epre) )
            f_log.write('%7.3f   %7.3f  %13.8f  %13.8f\n' % (alpha, hmax, energy, abs(energy-Econv))) 
            continue

        ''' Prepare an input file '''
        f_inp = open('tmp.inp', 'w')
        for  l in open(inp_tmpl):
            if l.find('cutoff_ele') != -1:
                f_inp.write('cutoff_ele = %10.5f\n' % max_cut)
            elif l.find('ewld_alpha') != -1:
                f_inp.write('ewld_alpha = %10.5f\n' % alpha)
            elif l.find('ewld_hmax') != -1:
                f_inp.write('ewld_hmax = %10.5f\n' % hmax)
            else:
                f_inp.write(l)
        f_inp.close()
    
        ''' Execute the code '''
        f_out = open('out','w')
        f_err = open('err','w')
        cmds = [ exe1, 'tmp.inp' ]
        subprocess32.call( cmds, stdout=f_out, stderr=f_err)
        f_out.close()
        f_err.close()
    
        ''' Read the Ewald energy output '''
        time = 0.0
        energy = 0.0
        flg_EWLD = True
        for l in open('out'):
            if flg_EWLD and l.find("EWLD: ") != -1:
                energy = float( l.strip().split()[-1] )
                flg_EWLD = False   ## To read only the first line of "EWLD" (corresponding the initial configuration)
            elif l.find('main_loop') != -1:
                time =  float( l.strip().split()[1] )
    
        data.append( (alpha, hmax, energy) )
        f_log.write('%7.3f   %7.3f  %13.8f  %13.8f   %9.5f\n' % (alpha, hmax, energy, abs(energy-Econv), time))
        f_log.flush()

        ''' For saving time
        if abs(energy - Econv) <= accuracy_for_exact:
            flg_skip = True
            Epre = Econv
            continue

        if abs(energy - Epre) <= accuracy_for_exact:
            flg_skip = True
        else:
            Epre = energy
        '''

    f_log.write('\n')

f_log.write('\n\n')
f_log.flush()

#data = [(0.05, 10.0, 0.0), (0.05,12.0, 0.0), (0.05,15.0,0.0)]


'''
For each (alpha, hmax),
First find the minimum cutoff_ele so that the energy-computation still achieves the accuracy.
and then measure the computation time using the minimum cutoff_ele.
'''
cut_series = []
i = 0 ; cut = max_cut
while (cut > 0.0):
    i += 1
    cut = max_cut - (5.0 * i)
    cut_series.append(cut)

fastest_time = 9999.9999
fastest_data = None

f_log.write('# alpha     Hmax    cutoff     T(main)    T(Neigh)   T(Energy)   T(Force)    T(E_Real)  T(E_Fourier) T(F_Real)  T(F_Fourier)\n')
for d in data:
    alpha  = d[0]
    hmax   = d[1]
    energy = d[2]

    if abs(energy - Econv) > accuracy:
        continue
    else:
        cutoff_ele = max_cut

    ''' Find the smallest possible cutoff_ele that still gives accurate energy '''
    for cut in cut_series:

        ''' Prepare an input file '''
        f_inp = open('tmp.inp', 'w')
        for  l in open(inp_tmpl):
            if l.find('cutoff_ele') != -1:
                f_inp.write('cutoff_ele = %10.5f\n' % cut)
            elif l.find('ewld_alpha') != -1:
                f_inp.write('ewld_alpha = %10.5f\n' % alpha)
            elif l.find('ewld_hmax') != -1:
                f_inp.write('ewld_hmax = %10.5f\n' % hmax)
            else:
                f_inp.write(l)
        f_inp.close()
    
        ''' Execute the code '''
        f_out = open('out','w')
        f_err = open('err','w')
        cmds = [ exe1, 'tmp.inp' ]
        subprocess32.call( cmds, stdout=f_out, stderr=f_err)
        f_out.close()
        f_err.close()
    
        ''' Read the Ewald energy output '''
        for l in open('out'):
            if l.find("EWLD: ") != -1:
                energy = float( l.strip().split()[-1] )
                break
    
        if abs(energy-Econv) > accuracy:
            break
        else:
            cutoff_ele = cut

    
    ''' Measure the computation time of short MD run (not single point energy) '''
    ''' Prepare an input file '''
    f_inp = open('tmp_md.inp', 'w')
    for  l in open(inp_tmpl_md):
        if l.find('cutoff_ele') != -1:
            f_inp.write('cutoff_ele = %10.5f\n' % cutoff_ele)
        elif l.find('ewld_alpha') != -1:
            f_inp.write('ewld_alpha = %10.5f\n' % alpha)
        elif l.find('ewld_hmax') != -1:
            f_inp.write('ewld_hmax = %10.5f\n' % hmax)
        else:
            f_inp.write(l)
    f_inp.close()
    
    ''' Execute the code '''
    f_out = open('out','w')
    f_err = open('err','w')
    cmds = [ exe2, 'tmp_md.inp' ]

    try:
        subprocess32.call( cmds, stdout=f_out, stderr=f_err, timeout=max_time)

    except subprocess32.TimeoutExpired as e:  # If it takes too long time
        f_out.close()
        f_err.close()
        tm_force = max_time
        tm_force_EwR = max_time
        tm_force_EwF = max_time
        tm_energy = max_time
        tm_energy_EwR = max_time
        tm_energy_EwF = max_time
        tm_main_loop = max_time

    else:
        f_out.close()
        f_err.close()

        ''' Read the timings output '''
        for l in open('out'):
            if l.find("_force(ele)") != -1:
                tm_force = float(l.split()[1])
            elif l.find("_force(ele_EwR)") != -1:
                tm_force_EwR = float(l.split()[1])
            elif l.find("_force(ele_EwF)") != -1:
                tm_force_EwF = float(l.split()[1])
            elif l.find("_neighbor(ele)") != -1:
                tm_neighbor = float(l.split()[1])
            elif l.find("_energy(ele)") != -1:
                tm_energy = float(l.split()[1])
            elif l.find("_energy(ele_EwR)") != -1:
                tm_energy_EwR = float(l.split()[1])
            elif l.find("_energy(ele_EwF)") != -1:
                tm_energy_EwF = float(l.split()[1])
            elif l.find("main_loop") != -1:
                tm_main_loop = float(l.split()[1])

    f_log.write('%7.3f   %7.3f   %5.1f     %9.5f   %9.5f  %9.5f  %9.5f    %9.5f  %9.5f    %9.5f  %9.5f\n' % 
           (alpha, hmax, cutoff_ele, 
            tm_main_loop,
            tm_neighbor,   tm_energy,  tm_force, 
            tm_energy_EwR, tm_energy_EwF, 
            tm_force_EwR,  tm_force_EwF ))
    f_log.flush()

    if tm_main_loop < fastest_time:
        fastest_time = tm_main_loop
        fastest_data = d

f_log.write('\n\n')
f_log.write('#Fastest time: %f\n' % fastest_time)
f_log.write('#Fastest parameters:\n')
f_log.write('#   alpha = %f\n' % fastest_data[0])
f_log.write('#   Hmax  = %f\n' % fastest_data[1])
f_log.close()
