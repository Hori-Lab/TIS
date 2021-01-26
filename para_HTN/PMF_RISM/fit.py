#!/usr/bin/env python

import sys
import numpy as np
from scipy.optimize import curve_fit
import math

Xs = []
Ys = []
for l in open('pmf_Mg_P_1264'):
    lsp = l.split()
    Xs.append(float(lsp[0]))
    Ys.append(float(lsp[1]))

#####################################################################
# Fit the short range
# r <= 2.25
#####################################################################
xx = []
yy = []
for x, y in zip(Xs,Ys):
    if x <= 2.25:
        xx.append(x)
        yy.append(y)

def repulsion(x, r0, beta):
    return np.power(r0 / x,  beta)

(r0, beta), dev = curve_fit(repulsion, xx, yy, p0=[2.0,12])

print ('#r0= %f' % (r0,))
print ('#beta= %f' % (beta,))
#a = 2.59972102
#b = 15.51478404


#####################################################################
# Fit the long range
# 12 < x
#####################################################################
xx = []
yy = []
for x, y in zip(Xs,Ys):
    if 12.0 <= x:
        xx.append(x)
        yy.append(y)
        
def yukawa(x, ay, ry):
    return - ay * np.exp(-ry*x)/x

(ay,ry), dev = curve_fit(yukawa, xx, yy)

print ('#aY= %f' % (ay,))
print ('#rY= %f' % (ry,))
#c = 8.48651135
#r = 0.01137363


#####################################################################
# Fit the middle
# 2 < x < 10
#####################################################################
xx = []
yy = []
for x, y in zip(Xs,Ys):
    if x < 2.25 or 10.0 < x:
        continue
    xx.append(x)
    yy.append(y)

#####################################################################
# power 4, 4 exponentials

def pmf_4_4exp(x, alpha1, x01, w1, 
                  alpha2, x02, w2,
                  alpha4, x04, w4, 
                  alpha5, x05, w5 ):
    return  ( repulsion(x, r0, beta) 
            - alpha1 * np.exp(-np.power(x-x01, 4.0) / w1)
            + alpha2 * np.exp(-np.power(x-x02, 2.0) / w2)
            + alpha4 * np.exp(-np.power(x-x04, 2.0) / w4)
            + alpha5 * np.exp(-np.power(x-x05, 2.0) / w5)
            + yukawa(x, ay, ry)
            )

p0 = [ 2.0, 2.2, 0.1,  # alpha1, x01, w1
       4.0, 3.1, 0.4,  # alpha2, x02, w2
       0.57, 6.0, 0.45, # alpha4, x04, w4
       0.1, 8.7, 0.7]     # alpha5, x05, w5
bounds = [ [0.001, 2.0, 0.001,
            0.001, 2.5, 0.001,
            0.001, 5.5, 0.001,
            0.001, 7.0, 0.001],
           [20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 7.0, 20.0,
            20.0, 9.5, 20.0]]

para_4_4exp, dev = curve_fit(pmf_4_4exp, xx, yy, p0=p0, bounds=bounds)

print('# 4_4exp')
print('1: %f %f %f' % (para_4_4exp[0], para_4_4exp[1], para_4_4exp[2]))
print('2: %f %f %f' % (para_4_4exp[3], para_4_4exp[4], para_4_4exp[5]))
print('4: %f %f %f' % (para_4_4exp[6], para_4_4exp[7], para_4_4exp[8]))
print('5: %f %f %f' % (para_4_4exp[9], para_4_4exp[10], para_4_4exp[11]))




#####################################################################
# power 4, 5 exponentials

def pmf_4_5exp(x, alpha1, x01, w1, 
                  alpha2, x02, w2,
                  alpha3, x03, w3,
                  alpha4, x04, w4, 
                  alpha5, x05, w5 ):
    return  ( repulsion(x, r0, beta) 
            - alpha1 * np.exp(-np.power(x-x01, 4.0) / w1)
            + alpha2 * np.exp(-np.power(x-x02, 2.0) / w2)
            + alpha3 * np.exp(-np.power(x-x03, 2.0) / w3)
            + alpha4 * np.exp(-np.power(x-x04, 2.0) / w4)
            + alpha5 * np.exp(-np.power(x-x05, 2.0) / w5)
            + yukawa(x, ay, ry)
            )

p0 = [ 2.0, 2.2, 0.1,  # alpha1, x01, w1
       4.0, 3.1, 0.4,  # alpha2, x02, w2
       2.2, 3.7, 0.3,  # alpha3, x03, w3
       0.57, 6.0, 0.45,# alpha4, x04, w4
       0.1, 8.7, 0.7]   # alpha5, x05, w5
bounds = [ [0.001, 2.0, 0.001,
            0.001, 2.5, 0.001,
            0.001, 3.2, 0.001,
            0.001, 5.5, 0.001,
            0.001, 7.0, 0.001],
           [20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 7.0, 20.0,
            20.0, 9.5, 20.0]]

para_4_5exp, dev = curve_fit(pmf_4_5exp, xx, yy, p0=p0, bounds=bounds)

print('# 4_5exp')
print('1: %f %f %f' % (para_4_5exp[0], para_4_5exp[1], para_4_5exp[2]))
print('2: %f %f %f' % (para_4_5exp[3], para_4_5exp[4], para_4_5exp[5]))
print('3: %f %f %f' % (para_4_5exp[6], para_4_5exp[7], para_4_5exp[8]))
print('4: %f %f %f' % (para_4_5exp[9], para_4_5exp[10], para_4_5exp[11]))
print('5: %f %f %f' % (para_4_5exp[12], para_4_5exp[13], para_4_5exp[14]))


#####################################################################
# power 4, 4 exponentials, with sigmoid

def pmf_4_4exp_sig(x, alpha1, x01, w1, 
                      alpha2, x02, w2,
                      alpha4, x04, w4, 
                      alpha5, x05, w5,
                      sig1, sig2):
    return  ( repulsion(x, r0, beta) / (1.0 + np.exp(sig1*(x-sig2)))
            - alpha1 * np.exp(-np.power(x-x01, 4.0) / w1)
            + alpha2 * np.exp(-np.power(x-x02, 2.0) / w2)
            + alpha4 * np.exp(-np.power(x-x04, 2.0) / w4)
            + alpha5 * np.exp(-np.power(x-x05, 2.0) / w5)
            + yukawa(x, ay, ry)
            )

p0 = [ 2.0, 2.2, 0.1,  # alpha1, x01, w1
       4.0, 3.1, 0.4,  # alpha2, x02, w2
       0.57, 6.0, 0.45,# alpha4, x04, w4
       0.1, 8.7, 0.7,  # alpha5, x05, w5
       1.0,1.0]
bounds = [ [0.001, 2.0, 0.001,
            0.001, 2.5, 0.001,
            0.001, 5.5, 0.001,
            0.001, 7.0, 0.001, 
            0.0, 0.0] ,
           [20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 7.0, 20.0,
            20.0, 9.5, 20.0,
            100.0, 3.0] ]

para_4_4exp_sig, dev = curve_fit(pmf_4_4exp_sig, xx, yy, p0=p0, bounds=bounds)

print('# 4_5exp_sig')
print('1: %f %f %f' % (para_4_4exp_sig[0], para_4_4exp_sig[1], para_4_4exp_sig[2]))
print('2: %f %f %f' % (para_4_4exp_sig[3], para_4_4exp_sig[4], para_4_4exp_sig[5]))
print('4: %f %f %f' % (para_4_4exp_sig[6], para_4_4exp_sig[7], para_4_4exp_sig[8]))
print('5: %f %f %f' % (para_4_4exp_sig[9], para_4_4exp_sig[10], para_4_4exp_sig[11]))
print('s: %f %f'    % (para_4_4exp_sig[12],para_4_4exp_sig[13]))

#####################################################################
# power 4, 5 exponentials, with sigmoid

def pmf_4_5exp_sig(x, alpha1, x01, w1, 
                      alpha2, x02, w2,
                      alpha3, x03, w3,
                      alpha4, x04, w4, 
                      alpha5, x05, w5,
                      sig1, sig2):
    return  ( repulsion(x, r0, beta) / (1.0 + np.exp(sig1*(x-sig2)))
            - alpha1 * np.exp(-np.power(x-x01, 4.0) / w1)
            + alpha2 * np.exp(-np.power(x-x02, 2.0) / w2)
            + alpha3 * np.exp(-np.power(x-x03, 2.0) / w3)
            + alpha4 * np.exp(-np.power(x-x04, 2.0) / w4)
            + alpha5 * np.exp(-np.power(x-x05, 2.0) / w5)
            + yukawa(x, ay, ry)
            )

p0 = [ 2.0, 2.2, 0.1,  # alpha1, x01, w1
       4.0, 3.1, 0.4,  # alpha2, x02, w2
       2.2, 3.7, 0.3,  # alpha3, x03, w3
       0.57, 6.0, 0.45,# alpha4, x04, w4
       0.1, 8.7, 0.7,  # alpha5, x05, w5
       1.0,1.0]
bounds = [ [0.001, 2.0, 0.001,
            0.001, 2.5, 0.001,
            0.001, 3.2, 0.001,
            0.001, 5.5, 0.001,
            0.001, 7.0, 0.001, 
            0.0, 0.0] ,
           [20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 7.0, 20.0,
            20.0, 9.5, 20.0,
            100.0, 3.0] ]

para_4_5exp_sig, dev = curve_fit(pmf_4_5exp_sig, xx, yy, p0=p0, bounds=bounds)

print('# 4_5exp_sig')
print('1: %f %f %f' % (para_4_5exp_sig[0], para_4_5exp_sig[1], para_4_5exp_sig[2]))
print('2: %f %f %f' % (para_4_5exp_sig[3], para_4_5exp_sig[4], para_4_5exp_sig[5]))
print('3: %f %f %f' % (para_4_5exp_sig[6], para_4_5exp_sig[7], para_4_5exp_sig[8]))
print('4: %f %f %f' % (para_4_5exp_sig[9], para_4_5exp_sig[10], para_4_5exp_sig[11]))
print('5: %f %f %f' % (para_4_5exp_sig[12], para_4_5exp_sig[13], para_4_5exp_sig[14]))
print('s: %f %f'    % (para_4_5exp_sig[15],para_4_5exp_sig[16]))

#####################################################################
# power 2, 4 exponentials

def pmf_2_4exp(x, alpha1, x01, w1, 
                  alpha2, x02, w2,
                  alpha4, x04, w4, 
                  alpha5, x05, w5 ):
    return  ( repulsion(x, r0, beta) 
            - alpha1 * np.exp(-np.power(x-x01, 2.0) / w1)
            + alpha2 * np.exp(-np.power(x-x02, 2.0) / w2)
            + alpha4 * np.exp(-np.power(x-x04, 2.0) / w4)
            + alpha5 * np.exp(-np.power(x-x05, 2.0) / w5)
            + yukawa(x, ay, ry)
            )

p0 = [ 2.0, 2.2, 0.1,  # alpha1, x01, w1
       4.0, 3.1, 0.4,  # alpha2, x02, w2
       0.57, 6.0, 0.45, # alpha4, x04, w4
       0.1, 8.7, 0.7]     # alpha5, x05, w5
bounds = [ [0.001, 2.0, 0.001,
            0.001, 2.5, 0.001,
            0.001, 5.5, 0.001,
            0.001, 7.0, 0.001],
           [20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 7.0, 20.0,
            20.0, 9.5, 20.0]]

para_2_4exp, dev = curve_fit(pmf_2_4exp, xx, yy, p0=p0, bounds=bounds)

print('# 2_4exp')
print('1: %f %f %f' % (para_2_4exp[0], para_2_4exp[1], para_2_4exp[2]))
print('2: %f %f %f' % (para_2_4exp[3], para_2_4exp[4], para_2_4exp[5]))
print('4: %f %f %f' % (para_2_4exp[6], para_2_4exp[7], para_2_4exp[8]))
print('5: %f %f %f' % (para_2_4exp[9], para_2_4exp[10], para_2_4exp[11]))


#####################################################################
# power 2, 5 exponentials

def pmf_2_5exp(x, alpha1, x01, w1, 
                  alpha2, x02, w2,
                  alpha3, x03, w3,
                  alpha4, x04, w4, 
                  alpha5, x05, w5 ):
    return  ( repulsion(x, r0, beta) 
            - alpha1 * np.exp(-np.power(x-x01, 2.0) / w1)
            + alpha2 * np.exp(-np.power(x-x02, 2.0) / w2)
            + alpha3 * np.exp(-np.power(x-x03, 2.0) / w3)
            + alpha4 * np.exp(-np.power(x-x04, 2.0) / w4)
            + alpha5 * np.exp(-np.power(x-x05, 2.0) / w5)
            + yukawa(x, ay, ry)
            )

p0 = [ 2.0, 2.2, 0.1,  # alpha1, x01, w1
       4.0, 3.1, 0.4,  # alpha2, x02, w2
       2.2, 3.7, 0.3,  # alpha3, x03, w3
       0.57, 6.0, 0.45,# alpha4, x04, w4
       0.1, 8.7, 0.7]   # alpha5, x05, w5
bounds = [ [0.001, 2.0, 0.001,
            0.001, 2.5, 0.001,
            0.001, 3.2, 0.001,
            0.001, 5.5, 0.001,
            0.001, 7.0, 0.001],
           [20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 7.0, 20.0,
            20.0, 9.5, 20.0]]

para_2_5exp, dev = curve_fit(pmf_2_5exp, xx, yy, p0=p0, bounds=bounds)

print('# 2_5exp')
print('1: %f %f %f' % (para_2_5exp[0], para_2_5exp[1], para_2_5exp[2]))
print('2: %f %f %f' % (para_2_5exp[3], para_2_5exp[4], para_2_5exp[5]))
print('3: %f %f %f' % (para_2_5exp[6], para_2_5exp[7], para_2_5exp[8]))
print('4: %f %f %f' % (para_2_5exp[9], para_2_5exp[10], para_2_5exp[11]))
print('5: %f %f %f' % (para_2_5exp[12], para_2_5exp[13], para_2_5exp[14]))


#####################################################################
# power 2, 4 exponentials, with sigmoid

def pmf_2_4exp_sig(x, alpha1, x01, w1, 
                      alpha2, x02, w2,
                      alpha4, x04, w4, 
                      alpha5, x05, w5,
                      sig1, sig2):
    return  ( repulsion(x, r0, beta) / (1.0 + np.exp(sig1*(x-sig2)))
            - alpha1 * np.exp(-np.power(x-x01, 2.0) / w1)
            + alpha2 * np.exp(-np.power(x-x02, 2.0) / w2)
            + alpha4 * np.exp(-np.power(x-x04, 2.0) / w4)
            + alpha5 * np.exp(-np.power(x-x05, 2.0) / w5)
            + yukawa(x, ay, ry)
            )

p0 = [ 2.0, 2.2, 0.1,  # alpha1, x01, w1
       4.0, 3.1, 0.4,  # alpha2, x02, w2
       0.57, 6.0, 0.45,# alpha4, x04, w4
       0.1, 8.7, 0.7,  # alpha5, x05, w5
       1.0,1.0]
bounds = [ [0.001, 2.0, 0.001,
            0.001, 2.5, 0.001,
            0.001, 5.5, 0.001,
            0.001, 7.0, 0.001, 
            0.0, 0.0] ,
           [20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 7.0, 20.0,
            20.0, 9.5, 20.0,
            100.0, 3.0] ]

para_2_4exp_sig, dev = curve_fit(pmf_2_4exp_sig, xx, yy, p0=p0, bounds=bounds)

print('# 2_4exp_sig')
print('1: %f %f %f' % (para_2_4exp_sig[0], para_2_4exp_sig[1], para_2_4exp_sig[2]))
print('2: %f %f %f' % (para_2_4exp_sig[3], para_2_4exp_sig[4], para_2_4exp_sig[5]))
print('4: %f %f %f' % (para_2_4exp_sig[6], para_2_4exp_sig[7], para_2_4exp_sig[8]))
print('5: %f %f %f' % (para_2_4exp_sig[9], para_2_4exp_sig[10], para_2_4exp_sig[11]))
print('s: %f %f'    % (para_2_4exp_sig[12],para_2_4exp_sig[13]))


#####################################################################
# power 2, 5 exponentials, with sigmoid

def pmf_2_5exp_sig(x, alpha1, x01, w1, 
                      alpha2, x02, w2,
                      alpha3, x03, w3,
                      alpha4, x04, w4, 
                      alpha5, x05, w5,
                      sig1, sig2):
    return  ( repulsion(x, r0, beta) / (1.0 + np.exp(sig1*(x-sig2)))
            - alpha1 * np.exp(-np.power(x-x01, 2.0) / w1)
            + alpha2 * np.exp(-np.power(x-x02, 2.0) / w2)
            + alpha3 * np.exp(-np.power(x-x03, 2.0) / w3)
            + alpha4 * np.exp(-np.power(x-x04, 2.0) / w4)
            + alpha5 * np.exp(-np.power(x-x05, 2.0) / w5)
            + yukawa(x, ay, ry)
            )

p0 = [ 2.0, 2.2, 0.1,  # alpha1, x01, w1
       4.0, 3.1, 0.4,  # alpha2, x02, w2
       2.2, 3.7, 0.3,  # alpha3, x03, w3
       0.57, 6.0, 0.45,# alpha4, x04, w4
       0.1, 8.7, 0.7,  # alpha5, x05, w5
       1.0,1.0]
bounds = [ [0.001, 2.0, 0.001,
            0.001, 2.5, 0.001,
            0.001, 3.0, 0.001,
            0.001, 5.5, 0.001,
            0.001, 6.0, 0.001, 
            0.0, 0.0] ,
           [20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 5.0, 20.0,
            20.0, 7.0, 20.0,
            20.0, 9.5, 20.0,
            100.0, 3.0] ]

para_2_5exp_sig, dev = curve_fit(pmf_2_5exp_sig, xx, yy, p0=p0, bounds=bounds)

print('# 2_5exp_sig')
print('1: %f %f %f' % (para_2_5exp_sig[0], para_2_5exp_sig[1], para_2_5exp_sig[2]))
print('2: %f %f %f' % (para_2_5exp_sig[3], para_2_5exp_sig[4], para_2_5exp_sig[5]))
print('3: %f %f %f' % (para_2_5exp_sig[6], para_2_5exp_sig[7], para_2_5exp_sig[8]))
print('4: %f %f %f' % (para_2_5exp_sig[9], para_2_5exp_sig[10], para_2_5exp_sig[11]))
print('5: %f %f %f' % (para_2_5exp_sig[12], para_2_5exp_sig[13], para_2_5exp_sig[14]))
print('s: %f %f'    % (para_2_5exp_sig[15],para_2_5exp_sig[16]))

print ('\n')


#####################################################################
#####################################################################
# For plot, and compute RMSE
#####################################################################

rmse_4_4exp = 0.0
rmse_4_4exp_sig = 0.0
rmse_4_5exp = 0.0
rmse_4_5exp_sig = 0.0
rmse_2_4exp = 0.0
rmse_2_4exp_sig = 0.0
rmse_2_5exp = 0.0
rmse_2_5exp_sig = 0.0
n = 0
for x, y in zip(Xs,Ys):
    y_4_4exp = pmf_4_4exp(x, *para_4_4exp)
    y_4_4exp_sig = pmf_4_4exp_sig(x, *para_4_4exp_sig)
    y_4_5exp = pmf_4_5exp(x, *para_4_5exp)
    y_4_5exp_sig = pmf_4_5exp_sig(x, *para_4_5exp_sig)
    y_2_4exp = pmf_2_4exp(x, *para_2_4exp)
    y_2_4exp_sig = pmf_2_4exp_sig(x, *para_2_4exp_sig)
    y_2_5exp = pmf_2_5exp(x, *para_2_5exp)
    y_2_5exp_sig = pmf_2_5exp_sig(x, *para_2_5exp_sig)
    print (x,y, 
           y_4_4exp    ,
           y_4_4exp_sig,
           y_4_5exp    ,
           y_4_5exp_sig,
           y_2_4exp    ,
           y_2_4exp_sig,
           y_2_5exp    ,
           y_2_5exp_sig)

    # Only compute errors in this range
    if x < 2.25 or 10.0 < x:
        continue
    n += 1
    rmse_4_4exp     += (y - y_4_4exp     )**2
    rmse_4_4exp_sig += (y - y_4_4exp_sig )**2
    rmse_4_5exp     += (y - y_4_5exp     )**2
    rmse_4_5exp_sig += (y - y_4_5exp_sig )**2
    rmse_2_4exp     += (y - y_2_4exp     )**2
    rmse_2_4exp_sig += (y - y_2_4exp_sig )**2
    rmse_2_5exp     += (y - y_2_5exp     )**2
    rmse_2_5exp_sig += (y - y_2_5exp_sig )**2

print ('\n')
print ('RMSE(4_4exp    ) = %f' % (math.sqrt(rmse_4_4exp/float(n)),))
print ('RMSE(4_4exp_sig) = %f' % (math.sqrt(rmse_4_4exp_sig/float(n)),))
print ('RMSE(4_5exp    ) = %f' % (math.sqrt(rmse_4_5exp/float(n)),))
print ('RMSE(4_5exp_sig) = %f' % (math.sqrt(rmse_4_5exp_sig/float(n)),))
print ('RMSE(2_4exp    ) = %f' % (math.sqrt(rmse_2_4exp/float(n)),))
print ('RMSE(2_4exp_sig) = %f' % (math.sqrt(rmse_2_4exp_sig/float(n)),))
print ('RMSE(2_5exp    ) = %f' % (math.sqrt(rmse_2_5exp/float(n)),))
print ('RMSE(2_5exp_sig) = %f' % (math.sqrt(rmse_2_5exp_sig/float(n)),))


##### Finatl: I chose 2_4exp_sig
f_out = open('fit_2_4exp_sig.out', 'w')
f_out.write('# r   PMF(r)   Fit(r)\n')
for x, y in zip(Xs,Ys):
    y_2_4exp_sig = pmf_2_4exp_sig(x, *para_2_4exp_sig)
    f_out.write('%f  %f  %f\n' % (x,y,y_2_4exp_sig,))
f_out.close()
