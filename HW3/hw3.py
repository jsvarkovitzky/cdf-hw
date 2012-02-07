#This program computes and plots the results of HW2 Question 1
#Jonathan Varkovitzky
#Jan 22, 2012

from numpy import *
from matplotlib import *
from pylab import *
import scipy as sp
import scipy.sparse
rc("text",usetex=True)

##################
## Load Airfoil ##
##################
def loadAirfoil():
    fname = 'body.dat'
    airfoilData = loadtxt(fname)
    return airfoilData

##################
## Plot Airfoil ##
##################
def plotAirfoil(airfoil):
    figure(1)
    plot(airfoil[:,1],airfoil[:,2])
    axis([0,1,-1,1])    
    title('NACA 0012 Airfoil')
    xlabel('x')
    ylabel('y')
    savefig('NACA0012')
    return

###################
## Initialize BC ##
###################
def initBC():
    x = zeros((n,m))
    y = zeros((n,m))
    
    #Read in points of airfoil as lower boundary
    for i in range(0,n):
        x[i,0] = airfoil[i,1]
        y[i,0] = airfoil[i,2]

    #Set outer BC to circle with a radius of 10 cord lenths
    dTheta = 2*pi/(n-1) #difference in angle between points
    for i in range(0,n):
        x[i,m-1] = r*cos(dTheta*i)
        y[i,m-1] = -r*sin(dTheta*i)
 
    #Set location of "cut"
    dx = (r - airfoil[m/2,1])/m #spacing from tail to boundary
    dx = 9./m
    for j in range(1,m-1):
        x[0,j] = 1+dx*j
        y[0,j] = 0
        x[n-1,j] = 1+dx*j
        y[n-1,j] = 0

    return (x,y)

 
################
## BC Plotter ##
################
def bcPlotter(x,y):
    figure(2)
    title('Mapping of Computaitonal Domain Boundaries to Physical Domain')
    xlabel('x')
    ylabel('y')
    plot(x[:,0],y[:,0],'r')
    plot(x[:,m-1],y[:,m-1],'g')
    plot(x[0,:],y[0,:],'b')
    plot(x[n-1,:],y[n-1,:],'k--')
    legend(('Bottom Boundary','Top Boundary','Left Boundary','Right Boundary'))
    
    return

##################
## Mesh Plotter ##
##################
def initMeshPlotter(x,y):
    figure(3)
    title('Mapping of Computaitonal Domain to Physical Domain')
    xlabel('x')
    ylabel('y')
    plot(x[:,:],y[:,:],'b.')
    plot(x[:,0],y[:,0],'r')
#    legend(('Transfinite Interpolation','NACA 0012 Airfoil'))
    savefig('initMesh')
    axis([-0.5,1.5,-0.5,0.5]) 
    savefig('initMeshZoom')
    return


#####################
## Initialize Mesh ##
#####################
def initMesh(x,y):
    Ax = zeros((n,m))
    Bx = zeros((n,m))
    Tx = zeros((n,m))
    Ay = zeros((n,m))
    By = zeros((n,m))
    Ty = zeros((n,m))
    for j in range(1,m):
        for i in range(1,n):
            L1j = float(j-m)/(0-m)
            L2j = float(j-0)/(m-0)
            L1i = float(i-n)/(0-n)
            L2i = float(i-0)/(n-0)
            Ax[i,j] = L1i*x[0,j]+L2i*x[n-1,j]
            Ay[i,j] = L1i*y[0,j]+L2i*y[n-1,j]
            Bx[i,j] = L1j*x[i,0]+L2j*x[i,m-1]
            By[i,j] = L1j*y[i,0]+L2j*y[i,m-1]
            Tx[i,j] = L1i*L1j*x[0,0]+L2i*L2j*x[n-1,m-1]+L1i*L2j*x[0,m-1]+L2i*L1j*x[n-1,0]
            Ty[i,j] = L1i*L1j*y[0,0]+L2i*L2j*y[n-1,m-1]+L1i*L2j*y[0,m-1]+L2i*L1j*y[n-1,0]

    Fx = Ax+Bx-Tx
    Fy = Ay+By-Ty
    #Insert Precomputed BCs
    Fx[:,0] = x[:,0]
    Fy[:,0] = y[:,0]
    Fx[:,m-1] = x[:,m-1]
    Fy[:,m-1] = y[:,m-1]
    Fx[0,:] = x[0,:]
    Fy[0,:] = y[0,:]
    Fx[n-1,:] = x[n-1,:]
    Fy[n-1,:] = y[n-1,:]
    return (Fx,Fy)
##################
## Main Program ##
##################

close ('all')
n = 129
m = 65
r = 10 #(10 cord lengths) x (cord lenght = 1) 
airfoil = loadAirfoil()
x,y = initBC()
#bcPlotter(x,y)
#plotAirfoil(airfoil)

x,y = initMesh(x,y)
initMeshPlotter(x,y)

show()
