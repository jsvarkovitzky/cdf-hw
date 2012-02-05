#This program computes and plots the results of HW2 Question 1
#Jonathan Varkovitzky
#Jan 22, 2012

from numpy import *
from matplotlib import *
from pylab import *
import scipy as sp
import scipy.sparse

#######################
## Analytic Solution ##
#######################

def analyticHeat(n,x,t):
    numTerms = 100
    uAnalytic = zeros(n+1)
    for k in range(1,numTerms+1):
#        if mod(k,2) == 1:
        uAnalytic = uAnalytic + 8*sin(k*pi/2)/(k*pi)**2*sin(k*pi*x)*exp(-(k*pi)**2*t)

    return uAnalytic


#############
## Set ICs ##
#############

def ICs(x,u):
    for i in range(0,n):
        
        if x[i]<0.5:
            u[i,0] = 2*x[i]
        else:
            u[i,0] = 2-2*x[i]
            
    return(u)

#############################
## Marches forward in time ##
#############################


def calcSoln(x,u,j):
    #Explicit Forward Difference Method
    if j == 0: 
        for k in range(0,T):
            for i in range(1,n):
                u[i,k+1] = u[i,k]+a*dt/dx**2*(u[i+1,k]-2*u[i,k]+u[i-1,k])


    #Implicit Backwards Difference Method
    if j == 1:
        #Define A matrix as stencil to act on internal points
        diag_rows = np.array([ones(n-1),-2*ones(n-1),ones(n-1)])
        positions = [-1, 0, 1]
        A = sp.sparse.spdiags(diag_rows, positions, n-1, n-1).todense()
        A = A*a*dt/dx**2
        w = u[1:n,0:T+1] #Selecting out only the interal points of u to act on

        #Update U over many time steps
        for k in range(0,T):
            w[:,k+1] = solve(eye(n-1)-A,w[:,k])
        u[1:n,0:T+1] = w[:,0:T+1] #Reincorperate internal points and BCs
    return(u)
##############
## Plotting ##
##############

def plotting(x,u,dt,i):
    uAnalytic = analyticHeat(n,x,i*dt)
    figure(i)
    clf()
    plot(x,u[:,i],'b')      #Plot computed solution
    plot(x,uAnalytic,'r--') #Plot analytic solution
    xlabel('x')
    ylabel('u(x,t)')
    legend(('Numerical Solution','Analytic Solution'))
    axis([0, 1, 0, 1])
    dT = int(dt*10000)
    if j == 0:
        savefig('heat_fdcd_%s_dt_%r'%(str(i).zfill(2),dT))
    elif j == 1:
        savefig('heat_bdcd_%s_dt_%r'%(str(i).zfill(2),dT))

##################
## Main Program ##
##################

close ('all')

# One step forward difference in time and centered difference in space
dx = 0.05
#dt = 0.0012
DT = [0.0012, 0.0013]
T = 50
xMin = 0
xMax = 1
n = int((xMax-xMin)/dx)

for m in range(0,2):
    dt = DT[m]
    print "Computing and plotting solutions for dt = %s" %dt
    for j in range(0,2):
        
        u = zeros((n+1,T+1))
        x = linspace(0,1,n+1)
        a = 1
        
        # Set ICs
        u = ICs(x,u)

        # Loop through remaining time steps to compute solution
        u = calcSoln(x,u,j)
        
        # Plot desired timestpes
        plotTimes = [0,1,10,50]
    
        for i in range(0,4):
            plotting(x,u,dt,plotTimes[i])



