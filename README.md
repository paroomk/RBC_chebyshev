This is the Chebyshev code.

To restart from a previous state saved in physical space, 
set the flag read_ICs = .true.  Make sure your restart files 
are saved as TICs.txt and VICs.txt for the temperature and 
vertical velocity fields respectively.  The data should be
saved such that the rows correspond to x-direction and 
columns correspond to y-direction.  Therefore, the j-th 
row and i-th column correspond to the field at y_j and 
x_i.

If you don't want to restart from a file the set 
read_ICs = .false. and specify your IC starting on line 
180.

Things to do:
1.) Input file
2.) VTK support?  .txt format may be enough 
3.) Adaptive time-stepping
