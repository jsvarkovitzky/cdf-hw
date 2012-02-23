from numpy import *
from scipy import *
from matplotlib import *
from pylab import *

################
## IC Plotter ##
################
def icPlotter(u,x):
    figure(1)
    title('rho IC')
    plot(x,u[0,:],'b')
    xlabel('x')
    ylabel('rho')
    
    figure(2)
    title('rho*u IC')
    plot(x,u[1,:],'b')
    xlabel('x')
    ylabel('rho*u')
    
    figure(3)
    title('rho*E IC')
    plot(x,u[2,:],'b')
    xlabel('x')
    ylabel('rho*E')
    
    show()
    
    return


###############
## u Plotter ##
###############
def uPlotter(u,x):
    figure(4)
    title('rho')
    plot(x,u[0,:],'b')
    xlabel('x')
    ylabel('rho')
    
    figure(5)
    title('rho*u')
    plot(x,u[1,:],'b')
    xlabel('x')
    ylabel('rho*u')
    
    figure(6)
    title('rho*E')
    plot(x,u[2,:],'b')
    xlabel('x')
    ylabel('rho*E')
    
    show()
    
    return

####################
## Predictor Step ##
####################
def uPredictor(u,f):
    uBar = zeros((3,n))
    uBar[:,1:n-1] = u[:,1:n-1]-tau*(f[:,2:n]-f[:,1:n-1])
    uBar[:,0] = u[:,0]
    uBar[:,-1] = u[:,-1]
    return uBar

####################
## Corrector Step ##
####################
def uCorrector(u,fBar):
    uBarBar = zeros((3,n))
    uBarBar[:,1:n-1] = u[:,1:n-1]-tau*(fBar[:,1:n-1]-fBar[:,0:n-2])
    uBarBar[:,0] = u[:,0]
    uBarBar[:,-1] = u[:,-1]
    return uBarBar

#################
## Update Step ##
#################
def uUpdate(u,uBar,uBarBar):
    u[:,:] = 1./2*(uBar[:,:]+uBarBar[:,:])
    return u
######################
## Compute F Vector ##
######################
def fFunc(u,f):
    p[:] = (g-1)*(u[2,:]-(u[1,:]**2/u[0,:])/2)#lec 12, pg 23
    f[0,:] = u[1,:]
    f[1,:] = u[1,:]**2/u[0,:]+p[:]
    f[2,:] = u[1,:]/u[0,:]*(u[2,:]+p[:])
    return f

##################
## Main Program ##
##################

close('all')

g = 1.4
n = 201
xmin = -0.5
xmax = 0.5
x = linspace(xmin,xmax,n)
dx = x[1]-x[0]
dt = dx/4000
tau = dt/dx

#I.C.s
pl = 1.*10**6
pr = 1.*10**5
p = zeros(n)
p[:] = where(x<0,pl,pr)
rhol = 1.0
rhor = 0.125

#Initialize u, uBar, f, fBar
u = zeros((3,n))
uBar = zeros((3,n))
f = zeros((3,n))
fBar = zeros((3,n))

#Fill u with I.C.s
u[0,:] = where(x<0,rhol,rhor)
u[1,:] = 0
u[2,:] = 1./(g-1)*p[:] #Based on Hirsch 1.4.1 & 1.4.25

#Plot I.C.s
#icPlotter(u,x)
N = 30
for i in range(0,N):
    f = fFunc(u,f)
    uBar = uPredictor(u,f)
    fBar = fFunc(uBar,fBar)
    uBarBar = uCorrector(u,fBar)
    u = uUpdate(u,uBar,uBarBar)


uPlotter(u,x)    


