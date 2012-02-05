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

#####################
## Initialize Mesh ##
#####################
def initMesh():
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
        y[i,m-1] = r*sin(dTheta*i)
 
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
    title('Mapping of Computaitonal Domain to Physical Domain')
    xlabel('x')
    ylabel('y')
    plot(x[:,0],y[:,0],'r')
    plot(x[:,m-1],y[:,m-1],'g')
    plot(x[0,:],y[0,:],'b')
    plot(x[n-1,:],y[n-1,:],'k--')
    legend(('Bottom Boundary','Top Boundary','Left Boundary','Right Boundary'))
    show()
    return

##################
## Main Program ##
##################

close ('all')
n = 129
m = 65
r = 10 #(10 cord lengths) x (cord lenght = 1) 
airfoil = loadAirfoil()
x,y = initMesh()
bcPlotter(x,y)
plotAirfoil(airfoil)





