# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 18:52:07 2015

@author: dsondak
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style

np.set_printoptions(threshold='inf')
matplotlib.style.use('fivethirtyeight')

data_BE = np.loadtxt('maxval_BE.txt')
data_IMEX = np.loadtxt('maxval_IMEX.txt')

time_IMEX = data_IMEX[:,0]
maxV_IMEX = data_IMEX[:,1]
maxT_IMEX = data_IMEX[:,2]

time_BE = data_BE[:,0]
maxV_BE = data_BE[:,1]
maxT_BE = data_BE[:,2]

fig1, ax1 = plt.subplots(1,1)
fig2, ax2 = plt.subplots(1,1)

alpha = 1.5585
gamma = np.pi/2.0
kappa = 0.1
decay_rate_T = kappa*(alpha**2.0 + gamma**2.0)
maxTe = maxT_BE[0]*np.exp(-decay_rate_T*time_BE)

gamma = 2.647427039234845
nu    = 0.05
decay_rate_V = nu*(alpha**2.0 + gamma**2.0)
maxVe = maxV_BE[0]*np.exp(-decay_rate_V*time_BE)

ax1.plot(time_BE,maxTe, label='analytical')
ax1.plot(time_IMEX,maxT_IMEX, ls = '--', label='IMEX-RK')
ax1.plot(time_BE,maxT_BE, ls = '--', label='Backward Euler')

ax1.legend()

fig1.savefig('heat-decay.pdf')

ax2.plot(time_BE,maxVe, label='analytical')
ax2.plot(time_IMEX,maxV_IMEX, ls = '--', label='IMEX-RK')
ax2.plot(time_BE,maxV_BE, ls = '--', label='Backward Euler')

ax2.legend()

fig2.savefig('stokes-decay.pdf')

plt.show()



