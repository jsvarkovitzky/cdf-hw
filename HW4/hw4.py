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
    uPred = zeros((3,n))
    uPred[:,1:n-1] = u[:,1:n-1] - tau*(f[:,2:n]-f[:,1:n-1])
    uPred[:,0] = u[:,0]
    uPred[:,n-1] = u[:,n-1]
    return uPred

####################
## Corrector Step ##
####################

#using method a for now
def uCorrect(u,uBar,fBar):
    uCor = zeros((3,n))
    uCor[:,1:n-1] = 1./2.*(u[:,1:n-1] + uBar[:,1:n-1]) - tau/2.*(fBar[:,1:n-1]-fBar[:,0:n-2])
    uCor[:,0] = u[:,0]
    uCor[:,n-1] = u[:,n-1]
    return uCor

##########################
## Compute the F vector ##
##########################

def fFunc(u,f):

    p = zeros(shape(u)[1])
    p[:] = (gamma-1)*(u[2,:]-u[1,:]**2/u[0,:]/2) #Lecture 12 pg 23
    f[0,:] = u[1,:]
    f[1,:] = u[1,:]**2/u[0,:] + p[:]
#    f[2,:] = u[1,:]*(gamma/(gamma-1)*p[:]/u[0,:]+(u[1,:]/u[0,:])**2/2)
#    f[2,:] = (u[2,:]+p[:])*u[1,:]/u[0,:]
#    f[2,:] = u[1,:]*u[2,:]/u[0,:] + p[:]/u[0,:]
    f[2,:] = u[1,:]/u[0,:]*(u[2,:]+p[:])
    return f

######################
## Newton Iteration ##
######################
def pNewton():
    R = 1
    P = 1 #Initial guess for P, keep small
    ul = 0
    ur = 0
    lim = 10**(-6)
    g = gamma
    a = alpha
    while R > lim:
        Pold = P
#        fp = sqrt(2/(g*(g-1)))*(P-1)/((1+a*P)**(1/2))-1/(g-1)*cl/cr*(1-(pr/pl*P)**((g-1)/(2*g)))+(ul-ur)/cr
#        fpPrime = sqrt(2/(gamma*(gamma-1)))*(1/(1+alpha*P)**(1/2)-(alpha*(P-1))/(1+alpha*P)**(3/2))+cl*pr/(cr*pl)/gamma*(pr/pl*P)**((gamma-1)/(2*gamma)-1)
        fp = (2/(g*(g-1)))**(1/2)*(P-1)/((1+a*P)**(1/2))-2/(g-1)*cl/cr*(1-(pr/pl*P)**((g-1)/(2*g)))-(ul-ur)/cr
        fpPrime = (2/(g*(g-1)))**(1/2)*(1/((1+a*P)**(1/2))-(a*(P-1))/(2*(1+a*P)**(3/2)))+cl/cr*pr/pl*1/g*(pr/pl*P)**((g-1)/(2*g)-1)
        P = P - fp/fpPrime
        R = abs(P-Pold)
#        print R

    return P

########################
## Ananlytic Solution ##
######################## 
def analyticSoln():
    P = pNewton()
    #Region 2
    p2 = P*pr
    u2 = (P-1)/(1+alpha*P)^(1/2)
    c2 = (P-1)*cr/(gamma*(u2-ur))+ur
    return
################
## Status Bar ##
################                                                                

def status(n,N):
    n = float(n)
    N = float(N)
    percent = n/N*100
    sys.stdout.write("[==> ]%3d%%\r" %percent)
    sys.stdout.flush()


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
pl = 1.*10**6
pr = 1.*10**5
rhol = 1.0
rhor = 0.125

alpha = (gamma+1)/(gamma-1)
cl = sqrt(gamma*pl/rhol)
cr = sqrt(gamma*pr/rhor)

dx = x[1]-x[0]
dt = dx/4/1000  ###Can set this to omptimize stability later###
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


figure(1)
title('Before')
plot(u[0,:],'b')
plot(u[1,:],'g')
plot(u[2,:]/10**7,'r')
show()



N = 30
for k in range(0,N):

    status(k,N)
    #calculate uBar
    uBar = uPredict(u,f)
    fBar = fFunc(uBar,fBar)
    u = uCorrect(u,uBar,fBar)
    f = fFunc(u,f)

figure(2)
title('After')
plot(u[0,:],'b')
plot(u[1,:],'g')
plot(u[2,:]/10**7,'r')
show()

P = pNewton()
print "p = %s" %P
