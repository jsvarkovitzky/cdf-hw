#This program solves the 1D Euler Equation for the shock tube problem
# via the MacCormick's method.

#Jonathan Varkovitzky
#2-13-2012








##################
## Main Program ##
##################

close('all')

figno = 0 #This is a counter to keep track of figure numbers

#Domain Properties
n = 201   #Number of grid points 
xmin = -0.5
xmax = 0.5

#I.C.s
pl = 10
pr = 1
rhol = 8
rhor = 1

