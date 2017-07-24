#!/bin/py

"""
"""

import numpy as np
import matplotlib.pyplot as plt
import config

plt.rcParams.update(config.pars)

def plot_ave_profile(f,y):
    """
    This function plots wall-normal profiles of fields.  The 
    field is first averaged over the x-direction.  Given a 
    field at a time t*, we plot \overline{f}(y,t*) where 
    \overline(.) denotes a horizontal average.

    The inputs to this function are:
       f:  The field to be plotted.  Rows are slices in y, 
           columns are slices in x
       y:  The wall-normal coordinate

    Note:  We should also allow the user to input a character 
           string indicating what field is being plotted (to 
           get the axes labels correct):.
    """
    favex = np.mean(f, axis=1) # Horizontal average

    # Plot
    fig, ax = plt.subplots(1,1, figsize=(15,9))
    ax.plot(favex, y)

    # Label
    ax.set_xlabel(r'$\overline{T}$')
    ax.set_ylabel(r'$y$')
    fig.tight_layout()

    # Save
    fig.savefig('Tave.pdf')
 
    plt.show()

if __name__ == "__main__":

   ### Load solution data ###
   Tyx = np.loadtxt('../Th/Tyx000015000.txt') # Temperature
   y = np.loadtxt('../y.txt') # y-coordinate
   Ny = np.shape(Tyx)[0] # Number of points in y direction
   Nx = np.shape(Tyx)[1] # Number of points in x direction

   y = np.hstack((1.0, y, -1.0)) # Build in top and bottom boundaries


   alpha = 1.5585 # Wavenumber we're working with
   x = np.linspace(-np.pi/alpha, np.pi/alpha, Nx+1) # Set up x-coordinate
                                                    # Note that we have added 
                                                    # 1 to account for the  
                                                    # repeated point (periodic bcs)
   
   X,Y = np.meshgrid(x,y) # Set up grid
   
   Tcond = -Y[:,0:-1] # Conduction solution
   
   # Top and bottom BCs for solution
   Ttop = np.zeros(Nx)
   Tbot = np.zeros(Nx)

   # Build top and bottom BCs into temperature field
   Tyx = np.vstack((Ttop, Tyx, Tbot))
   
   # Add conduction solution back in to field
   Tyx = Tyx + Tcond

   # Plot the horizontally-averaged temperature profile
   plot_ave_profile(Tyx, y)
