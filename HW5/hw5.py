#This program solves the 2D Euler Equations about the NACA-0012 airfoil
#using the Jameson Finite Volume Method

#Jonathan Varkovitzky
#Feb 29, 2012

from numpy import *
from scipy import *
from pylab import *
import time
import os


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
    sxx[0:n,0:m] = -ys[1:n+1,1:m+1]-ys[1:n+1,0:m]#/sqrt((ys[1:n+1,1:m+1]-ys[1:n+1,0:m])**2+(xs[1:n+1,1:m+1]-xs[1:n+1,0:m])**2)
    sxy[0:n,0:m] = (xs[1:n+1,1:m+1]-xs[1:n+1,0:m])#/sqrt((ys[1:n+1,1:m+1]-ys[1:n+1,0:m])**2+(xs[1:n+1,1:m+1]-xs[1:n+1,0:m])**2)
    syx[0:n,0:m] = sxy[0:n,0:m]
    syy[0:n,0:m] = -sxx[0:n,0:m]


    return (sxx,sxy,syx,syy)

################
## F Function ##
################
def f_func(u):
    F = zeros((4,n,m))
    p = zeros((n,m))

    p[:,:] = (g-1)*(u[3,:,:]-(u[1,:,:]**2+u[2,:,:]**2)/(2*u[0,:,:]))
    F[0,:,:] = u[1,:,:]
    F[1,:,:] = u[1,:,:]**2/u[0,:,:] + p[:,:]
    F[2,:,:] = u[1,:,:]*u[2,:,:]/u[0,:,:]
    F[3,:,:] = (u[3,:,:]+p[:,:])*(u[1,:,:]/u[0,:,:])

    return (F)

################
## G Function ##
################
def g_func(u):
    G = zeros((4,n,m))
    p = zeros((n,m))

    p[:,:] = (g-1)*(u[3,:,:]-(u[1,:,:]**2+u[2,:,:]**2)/(2*u[0,:,:]))
    G[0,:,:] = u[1,:,:]
    G[1,:,:] = u[1,:,:]*u[2,:,:]/u[0,:,:]
    G[2,:,:] = u[2,:,:]**2/u[0,:,:] + p[:,:]
    G[3,:,:] = (u[3,:,:]+p[:,:])*(u[2,:,:]/u[0,:,:])

    return (G)
######################
## Flux Calculation ##
######################
def flux(u,sx,sy):

    

    return Flux

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

###########################
## Normal Vector Plotter ##
###########################
def normalPlotter(x,y,xs,ys,sxx,sxy,syx,syy):
    figure(2)
    for i in range(0,n+1):
        plot(xs[i,:],ys[i,:],'b')
    for j in range(0,m+1):
        plot(xs[:,j],ys[:,j],'b')
    scale = 10000
 #   plot(x,y,'r.')
    quiver(x,y,sxx*scale,sxy*scale)
    quiver(x,y,syx*scale,syy*scale)
    title('Cell Normal Directions')
    xlabel('x')
    ylabel('y')

    return
#############################
## Velocity Vector Plotter ##
#############################
def velocityPlotter(x,y,xs,ys,x_vec,y_vec):
    figure(3)
    for i in range(0,n+1):
        plot(xs[i,:],ys[i,:],'b')
    for j in range(0,m+1):
        plot(xs[:,j],ys[:,j],'b')
 #   plot(x,y,'r.')
    quiver(x,y,x_vec,y_vec)
    title('Velocity Directions')
    xlabel('x')
    ylabel('y')

    return

##################
## Main Program ##
##################
start_time = time.time()
n = 128
m = 64
g = 1.4
p_0 = 10**5
rho = 1
M_stream = 0.85

print "Loading Grid Points..."
(x_mesh,y_mesh) = loadXY()
print "Computing Cell Centers..."
(x,y) = cellCenters(x_mesh,y_mesh)
print "Computing Cell Areas..."
omega = cellAreas(x_mesh,y_mesh)
print "Computing Cell Normals..."
(sxx,sxy,syx,syy) = cellNorms(x_mesh,y_mesh)
#meshPlotter(x_mesh,y_mesh,x,y)
#normalPlotter(x,y,x_mesh,y_mesh,sxx,sxy,syx,syy)

print "Initializing vectors..."
#Set IC's and BC's together assuming an initial uniform velocity field
u = zeros((4,n,m))
F = zeros((4,n,m))
G = zeros((4,n,m))
u[0,:,:] = 1.0*1000 #initialize rho
u[1,:,:] = M_stream#initialize x velocity
u[2,:,:] = 0#initialize y velocity
u[3,:,:] = p_0/(g-1)+rho*(u[1,:,:]**2+u[2,:,:]**2)/2#initialize energy
#Set initial F and G values
F = f_func(u)
G = g_func(u) 
#Set initial Boundaries

#Branch Cut
u[:,-1,:] = u[:,1,:]
u[:,0,:] = u[:,-2,:]

#Airfoil
u[2,:,0] = -u[2,:,0]#negate y velocity



print time.time() - start_time, "seconds"
