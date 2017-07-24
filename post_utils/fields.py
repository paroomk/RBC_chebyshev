#!/bin/py

import numpy as np
import matplotlib.pyplot as plt

### Load solution data ###
Tyx = np.loadtxt('Th/Tyx000007000.txt') # Temperature
Vyx = np.loadtxt('VV/Vyx000007000.txt') # Velocity
y = np.loadtxt('y.txt') # y-coordinate
Ny = np.shape(Tyx)[0] # Number of points in y direction
Nx = np.shape(Tyx)[1] # Number of points in x direction

print Nx, Ny

y = np.hstack((1.0, y, -1.0)) # Build in top and bottom boundaries

alpha = 1.5585 # Wavenumber we're working with
x = np.linspace(-np.pi/alpha, np.pi/alpha, Nx+1) # Set up x-coordinate
                                                 # Note that we have added 
                                                 # 1 to account for the  
                                                 # repeated point (periodic bcs)

X,Y = np.meshgrid(x,y) # Set up grid

Tcond = -Y # Conduction solution

# Top and bottom BCs for solution
Ttop = np.zeros(Nx)
Tbot = np.zeros(Nx)
Vtop = np.zeros(Nx)
Vbot = np.zeros(Nx)

# Build top and bottom BCs into temperature field
Tyx = np.vstack((Ttop, Tyx, Tbot))
Tyx = np.c_[Tyx, Tyx[:,0]]

# Add conduction solution back in to field
Tyx = Tyx + Tcond

# Plot temperature contours
figT, axT = plt.subplots(1,1)
axT.contourf(X, Y, Tyx)

# Build top and bottom BCs into velocity field
Vyx = np.vstack((Vtop, Vyx, Vbot))
Vyx = np.c_[Vyx, Vyx[:,0]]

# Plot vertical velocity contours
figV, axV = plt.subplots(1,1)
axV.contourf(X, Y, Vyx)

plt.show()
