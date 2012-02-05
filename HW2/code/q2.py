#This program computes and plots the results of HW2 Question 2
#Jonathan Varkovitzky
#Jan 25, 2012

from numpy import *
from matplotlib import *
from pylab import *
import scipy as sp
import scipy.sparse

#######################
## Analytic Solution ##
#######################

def advectAnalytic(n,x,u,i,dt):
    t = i*dt
    uAnalytic = zeros(n+1)

    for i in range(0,n):
        if (x[i] < -5+a*t) and (x[i]>-7+a*t):
            uAnalytic[i] = 1

    u[:,dt] = uAnalytic

    return u

#############################
## Marches forward in time ##
#############################

def calcSoln(x,u,j):
    sigma = a*dt/dx
#    print "sigma = %s" %sigma

    #Explicit Forward Difference Method

    if j == 0: 
        for k in range(0,T):
            for i in range(1,n):
                u[i,k+1] = u[i,k]-sigma*(u[i,k]-u[i-1,k])


    #Lax-Wendroff Method
    if j == 1:
        for k in range(0,T):
            for i in range(1,n):
                u[i,k+1] = u[i,k]-sigma/2*(u[i+1,k]-u[i-1,k])+sigma**2/2*(u[i+1,k]-2*u[i,k]+u[i-1,k])
        
    #Implicit Method of Choice
    if j == 2:
        #Define A matrix as stencil to act on internal points
        diag_rows = np.array([sigma/2*ones(n-1),zeros(n-1),-sigma/2*ones(n-1)])
        positions = [-1, 0,1]
        A = sp.sparse.spdiags(diag_rows, positions, n-1, n-1).todense()
        w = u[1:n,0:T+1] #Selecting out only the interal points of u to act on

        #Update U over many time steps
        for k in range(0,T):
            w[:,k+1] = solve(eye(n-1)-A,w[:,k])
        u[1:n,0:T+1] = w[:,0:T+1] #Reincorperate internal points and BCs

    return(u)
#############
## Errrors ##
#############
def errorCalc(sigma,cfl,j):
    phi = linspace(0.001,pi,100)
    if j == 0:
        G = abs(1-sigma*(1-exp(-1j*phi)))
        PHI = arctan((sigma*phi)/(1-sigma+sigma*cos(phi)))/(sigma*phi)
        meth = 'fdcd'
    elif j == 1:
        G = abs(1-1j*sigma*sin(phi)+sigma**2*cos(phi)-sigma**2)
        PHI = arctan((sigma*phi)/(1-sigma**2*(cos(phi)-1)))/(sigma*phi)
        meth = 'LW'
    elif j == 2:
        G = abs(1-sigma*1j*sin(phi))
        PHI = arctan(sigma*sin(phi))/(sigma*phi)
        meth = 'Implicit'

    print "error for j = %s"%j
    figure(10)
    clf()
    plot(phi,G)
    cflNum = int(cfl*10*2)
    xlabel('Phase')
    ylabel('Dispersion Error')
    savefig('diff_error_cfl_%r_%s' %(cflNum,meth))

    figure(11)
    clf()
    plot(phi,PHI)
    xlabel('Phase')
    ylabel('Dissipation Error')
    savefig('disp_error_cfl_%r_%s' %(cflNum,meth))              
                 
##############
## Plotting ##
##############

def plotting(x,u,dt,cfl,i):
    uAnalytic = advectAnalytic(n,x,u,i,dt)
    figure(i)
    clf()
    plot(x,u[:,i],'b')      #Plot computed solution
    plot(x,uAnalytic[:,dt],'r--') #Plot analytic solution
    xlabel('x')
    ylabel('u(x,t)')
    legend(('Numerical Solution','Analytic Solution'))
    axis([-10, 10, -0.3, 1.3])
    dT = int(dt*10000)
    if j == 0:
        savefig('advect_fdcd_%s_dt_%r'%(str(i).zfill(2),dT))
    elif j == 1:
        savefig('advect_lw_%s_dt_%r'%(str(i).zfill(2),dT))
    elif j == 2:
        savefig('advect_implicit_%s_dt_%r'%(str(i).zfill(2),dT))
##################
## Main Program ##
##################

close ('all')

# One step forward difference in time and centered difference in space
dx = 0.05
a = 0.5
#dt = 0.0012
CFL = array([0.6, 1.2])
DT = CFL*dx/a
xMin = -10
xMax = 10
n = int((xMax-xMin)/dx)

for m in range(0,2):
    dt = DT[m]
    T = int(15/dt)
    cfl = CFL[m]

    print "Computing and plotting solutions for CFL = %s" %cfl
    for j in range(0,3):

        u = zeros((n+1,T+1))
        x = linspace(xMin,xMax,n+1)
        # Set ICs
        u = advectAnalytic(n,x,u,0,0)
        # Loop through remaining time steps to compute solution
        u = calcSoln(x,u,j)

        # Plot desired timestpes
        plotTimes = [T]
    
        for i in range(0,1):
            plotting(x,u,dt,cfl,plotTimes[i])
        
        sigma = a*dt/dx
        errorCalc(sigma/2,cfl/2,j)    

