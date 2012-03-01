#This program solves the 2D Euler Equations about the NACA-0012 airfoil
#using the Jameson Finite Volume Method

#Jonathan Varkovitzky
#Feb 29, 2012

from numpy import *
from scipy import *
from pylab import *
################
## Load Grids ##
################
def loadXY():
    x = loadtxt('x_pts.txt')
    y = loadtxt('y_pts.txt')
    return (x,y)

##########################
## Compute Cell Centers ##
##########################
def cellCenters(xs,ys):
    x = zeros((n,m))
    y = zeros((n,m))
    for k in range(0,m):
        x[:,k] = 1./4*(xs[0:n,k]+xs[1:n+1,k]+xs[0:n,k+1]+xs[1:n+1,k+1])
        y[:,k] = 1./4*(ys[0:n,k]+ys[1:n+1,k]+ys[0:n,k+1]+ys[1:n+1,k+1])
    return (x,y)

########################
## Compute Cell Areas ##
########################
def cellAreas(xs,ys):
    omega = zeros((n,m))
    for i in range(0,n):
        for j in range(0,m):
            xac = array([xs[i,j+1]-xs[i+1,j],ys[i,j+1]-ys[i+1,j]])
            xbd = array([xs[i+1,j+1]-xs[i,j],ys[i+1,j+1]-ys[i,j]])
            omega[i,j] = 1./2*abs(cross(xac,xbd))
    return omega

##########################
## Compute Cell Normals ##
##########################
def cellNorms(xs,ys):
    sxx = zeros((n,m))
    sxy = zeros((n,m))
    syx = zeros((n,m))
    syy = zeros((n,m))
    
    sxx[0:n,0:m] = ys[1:n+1,1:m+1]-ys[1:n+1,0:m]
    sxy[0:n,0:m] = -(xs[1:n+1,1:m+1]+xs[1:n+1,0:m])
    syx[0:n,0:m] = sxy[0:n,0:m]
    syy[0:n,0:m] = -sxx[0:n,0:m]

    return (sxx,sxy,syx,syy)

##################
## Mesh Plotter ##
##################
def meshPlotter(xs,ys,x,y):
    figure(1)
    for i in range(0,n+1):
        plot(xs[i,:],ys[i,:],'b')
    for j in range(0,m+1):
        plot(xs[:,j],ys[:,j],'b')
    plot(x,y,'r.')
    xlabel('x')
    ylabel('y')
    title('Plot with Mesh and Cell Centers')
    return




##################
## Main Program ##
##################

n = 129-1
m = 65-1

print "Loading Grid Points..."
(x_mesh,y_mesh) = loadXY()
print "Computing Cell Centers..."
(x,y) = cellCenters(x_mesh,y_mesh)
print "Computing Cell Areas..."
omega = cellAreas(x_mesh,y_mesh)
print "Computing Cell Normals..."
(sxx,sxy,syx,syy) = cellNorms(x_mesh,y_mesh)
#meshPlotter(x_mesh,y_mesh,x,y)
#show()
