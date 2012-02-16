#This program solves the 1D Euler Equation for the shock tube problem
# via the MacCormick's method.

#Jonathan Varkovitzky
#2-13-2012

from numpy import *
from scipy import *
from matplotlib import *
from pylab import *

####################
## Predictor Step ##
####################

#using method a for now
def uPredict(u,f):
    uBar[:,1:n-2] = u[:,1:n-2] - tau*(f[:,2:n-1]-f[:,1:n-2])
    uBar[:,0] = u[:,0]
    uBar[:,n-1] = u[:,n-1]
    print "uBar:"
    print uBar

    return uBar

####################
## Corrector Step ##
####################

#using method a for now
def uCorrect(u,uBar,fBar):
    uBar[:,1:n-1] = 1./2*(u[:,1:n-1] + uBar[:,1:n-1]) - tau/2*(fBar[:,1:n-1]-fBar[:,0:n-2])
    return uBar

##########################
## Compute the F vector ##
##########################

def fFunc(u,f):
    
    p = (u[2,:]-u[1,:]**2/2)*(gamma-1)/u[0,:]
    f[0,:] = u[1,:]
    f[1,:] = u[1,:]**2/u[0,:] + p[:]
    f[2,:] = u[1,:]*(gamma/(gamma-1)*p[:]/u[0,:]+(u[1,:]/u[0,:])**2/2)
    return f

##################
## Main Program ##
##################

close('all')

#Various Constants
figno = 0   #This is a counter to keep track of figure numbers
gamma = 1.4
R = 8.3145
Tl = 273
Tr = Tl*0.8
kEnd = 100

#Domain Properties
n = 201   #Number of grid points 
xmin = -0.5
xmax = 0.5
x = linspace(xmin,xmax,n)

#I.C.s
pl = 10.0*10**5
pr = 1.0*10**5
rhol = 8.
rhor = 1.

dx = x[1]-x[0]
dt = dx/4  ###Can set this to omptimize stability later###
tau = (dt/dx)
#Initialize the vectors for storing data
#Note a history is not needed for the assignment and therefore is not kept
u = zeros((3,n))
uBar = zeros((3,n))
f = zeros((3,n))
fBar = zeros((3,n))
uNew = zeros((3,n))
fNew = zeros((3,n))


#Initiate I.C.s
u[0,:] = where(x<=0, rhol,rhor)
u[1,:] = 0.0
u[2,:] = where(x<0,1/(gamma-1)*u[0,:]*pl+u[1,:]**2/2,1/(gamma-1)*u[0,:]*pr+u[1,:]**2/2)

f = fFunc(u,f)

for k in range(0,10):

    #calculate uBar
    uBar = uPredict(u,f)
    fBar = fFunc(uBar,f)
    u = uCorrect(u,uBar,fBar)
    f = fFunc(u,f)

