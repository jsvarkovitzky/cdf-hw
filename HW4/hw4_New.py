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
######################
## Newton Iteration ##
######################
def pNewton():
    R = 1
    P = 1 #Initial guess for P, keep small
    lim = 10**(-6)
    while R > lim:
        Pold = P
        fp = (2./(g*(g-1)))**(1/2)*(P-1)/(1+a*P)**(1/2)-2/(g-1)*cl/cr*(1-(pr/pl*P)**((g-1)/2*g))+(ul-ur)/cr
        fpPrime = (2./(g*(g-1)))**(1/2)*(1/(1+a*P)**(1/2)-(a*(P-1))/(2*(1+a*P)**(3/2)))+cl/cr*pr/pl*(pr/pl*P)**((-g-1)/(2*g))
        P = P - fp/fpPrime
        R = abs(P-Pold)
    return P

###############################
## Compute Analytic Solution ##
###############################
def computeAnalytic(ts):
    
    p = zeros(n)
    u = zeros(n)
    rho = zeros(n)
    C = zeros(n)
    H = zeros(n)
    M = zeros(n)

    #region 5 (note no region 4)
    u[:] = 2/(g+1)*(x[:]/ts+cl+(g-1)/2*ul)
    C[:] = u[:]-x[:]/ts
    p[:] = pl*(C[:]/cl)**(2*g/(g-1))
    rho[:] = rhol*(p[:]/pl)**(1/g)
    M[:] = u[:]/C[:]
    H[:] = g/(g-1)*p[:]/rho[:]+u[:]**2/2

    #region 2
    p2 = P*pr
    u2 = (P-1)/sqrt(1+a*P)*cr/sqrt(g*(g-1)/2)
    rho2 = (1+a*P)/(a+P)*rhor
    C2 = (P-1)*cr/(g*(u2-ur))+ur
    H2 = g/(g-1)*p2/rho2+u2**2/2
    M2 = u2/C2

    #region 3
    p3 = p2  
    u3 = u2
    rho3 = rhol*(p3/pl)**(1/g)
    C3 = sqrt(g*p3/rho3)
    H3 = g/(g-1)*p3/rho3+u3**2/2
    M3 = u3/C3
    v = u2
   
    #fill values into vectors of region 3
    u[:] = where(x>v*ts,u3,u)
    p[:] = where(x>v*ts,p3,p)
    rho[:] = where(x>v*ts,rho3,rho)
    H[:] = where(x>v*ts,H3,H)
    C[:] = where(x>v*ts,C3,C)
    M[:] = where(x>v*ts,M3,M)
    
    #fill values into vectors of region 2
    CS = (P-1)*cr**2/(g*(u2-ur))+ur

    u[:] = where(x>CS*ts,u2,u)
    p[:] = where(x>CS*ts,p2,p)
    rho[:] = where(x>CS*ts,rho2,rho)
    H[:] = where(x>CS*ts,H2,H)
    C[:] = where(x>CS*ts,C2,C)
    M[:] = where(x>CS*ts,M2,M)
    
    #fill in values in left region
    left_lim = -ts*((g-1)/2*ul+cl)
    u[:] = where(x<left_lim,ul,u)
    p[:] = where(x<left_lim,pl,p)
    rho[:] = where(x<left_lim,rhol,rho)
    Hl = g/(g-1)*pl/rhol+ul**2/2
    Ml = ul/cl               
    H[:] = where(x<left_lim,Hl,H)
    C[:] = where(x<left_lim,cl,C)
    M[:] = where(x<left_lim,Ml,M)
    
    figure(17)
    plot(x,rho)

    print CS*ts, v*ts

    return

##################
## Main Program ##
##################

close('all')

g = 1.4
a = (g+1)/(g-1)
n = 201
xmin = -0.5
xmax = 0.5
x = linspace(xmin,xmax,n)
dx = x[1]-x[0]
#dt = dx/40000
#tau = dt/dx

#I.C.s
pl = 1.*10**5
pr = 1.*10**4
p = zeros(n)
p[:] = where(x<0,pl,pr)
rhol = 1.0
rhor = 0.125
ul = 0
ur = 0

#Initialize sound speeds
cl = sqrt(g*pl/rhol)
cr = sqrt(g*pr/rhor)

#Initialize u, uBar, f, fBar
u = zeros((3,n))
uBar = zeros((3,n))
f = zeros((3,n))
fBar = zeros((3,n))

#Fill u with I.C.s
u[0,:] = where(x<0,rhol,rhor)
u[1,:] = 0
u[2,:] = 1./(g-1)*p[:] #Based on Hirsch 1.4.1 & 1.4.25

rho = zeros(n)
CS = zeros(n)
#Plot I.C.s
#icPlotter(u,x)
N = 100
for i in range(0,N):
    p[:] = (g-1)*(u[2,:]-(u[1,:]**2/u[0,:])/2)
    rho[:] = u[0,:]
    C = sqrt(g*p[:]/rho[:])
    CS = u[1,:]/u[0,:] + C[:]
    amax = max(CS)
    dt = 0.95*dx/amax
    tau = dt/dx

    f = fFunc(u,f)
    uBar = uPredictor(u,f)
    fBar = fFunc(uBar,fBar)
    uBarBar = uCorrector(u,fBar)
    u = uUpdate(u,uBar,uBarBar)


#uPlotter(u,x)    

#P = pNewton()
P = 3.0312995
ts = 6.7479*10**(-4)
computeAnalytic(ts)

