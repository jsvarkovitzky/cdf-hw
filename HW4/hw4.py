#This program solves the 1D Euler Equation for the shock tube problem
# via the MacCormick's method.

#Jonathan Varkovitzky
#2-13-2012

from numpy import *
from scipy import *
from matplotlib import *
from pylab import *

##################
## Main Program ##
##################

close('all')

figno = 0 #This is a counter to keep track of figure numbers

#Domain Properties
n = 201   #Number of grid points 
xmin = -0.5
xmax = 0.5
x = linspace(xmin,xmax,n)

#I.C.s
pl = 10
pr = 1
rhol = 8
rhor = 1

#Initialize the vector for u to store p and rho in rows 0 and 1 respectively
#******POSSIBLE THE INCORRECT FORM OF U AND NEED TO HAVE AN F AS WELL****
u = zeros((2,n))
u[0,:] = where(x<=0, pl,u[0,:])
u[0,:] = where(x>0,  pr,u[0,:])
u[1,:] = where(x<=0, rhol,u[1,:])
u[1,:] = where(x>0,  rhor,u[1,:])

