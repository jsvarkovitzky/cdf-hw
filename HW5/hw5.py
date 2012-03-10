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
    #cell centers will have one less point than the grid
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
    omega = zeros((N,M))
    for i in range(0,N-1):
        for j in range(0,M-1):
            xac = array([xs[i,j+1]-xs[i+1,j],ys[i,j+1]-ys[i+1,j]])
            xbd = array([xs[i+1,j+1]-xs[i,j],ys[i+1,j+1]-ys[i,j]])
            omega[i,j] = 1./2*abs(cross(xac,xbd))
    return omega

##########################
## Compute Cell Normals ##
##########################
def cellNorms(xs,ys):
    #Note that radial normals are going clockwise, this is correct
    #as our loops go through idicies clockwise and require the normals
    #to be in the same direction.  To test directionality remove indicies
    #in the mesh plotting as desired.

    sxx = zeros((n,m))
    sxy = zeros((n,m))
    syx = zeros((n,m))
    syy = zeros((n,m))
    sxx[0:n,0:m] = ys[1:N,1:M]-ys[1:N,0:M-1]
    sxy[0:n,0:m] = -(xs[1:N,1:M]-xs[1:N,0:M-1])
    syx[0:n,0:m] = sxy[:,:]
    syy[0:n,0:m] = -sxx[:,:]


    return (sxx,sxy,syx,syy)

################
## F Function ##
################
def f_func(u):
    N = shape(u)[1]
    M = shape(u)[2]
    F = zeros((4,N,M))
    p = zeros((N,M))
    
    u[0,:,:] = rho

    p[:,:] = (g-1)*(u[3,:,:]-(u[1,:,:]**2+u[2,:,:]**2)/(2*u[0,:,:]))
    F[0,:,:] = u[1,:,:]
    F[1,:,:] = u[1,:,:]**2/u[0,:,:] + p[:,:]
    F[2,:,:] = u[1,:,:]*u[2,:,:]/u[0,:,:]
    F[3,:,:] = (u[3,:,:]+p[:,:])*(u[1,:,:]/u[0,:,:])

    #Impose airfoil BCs
    F[0,:,0] = 0
    F[1,1:-1,0] = p[1:-1,0]*sxx[:,0]
    F[2,1:-1,0] = p[1:-1,0]*sxy[:,0]
    F[3,:,0] = 0
    return (F)

################
## G Function ##
################
def g_func(u):

    N = shape(u)[1]
    M = shape(u)[2]
    G = zeros((4,N,M))
    p = zeros((N,M))

    u[0,:,:] = rho

    p[:,:] = (g-1)*(u[3,:,:]-(u[1,:,:]**2+u[2,:,:]**2)/(2*u[0,:,:]))
    G[0,:,:] = u[1,:,:]
    G[1,:,:] = u[1,:,:]*u[2,:,:]/u[0,:,:]
    G[2,:,:] = u[2,:,:]**2/u[0,:,:] + p[:,:]
    G[3,:,:] = (u[3,:,:]+p[:,:])*(u[2,:,:]/u[0,:,:])

    #Impose airfoil BCs
    G[0,:,0] = 0
    G[1,1:-1,0] = p[1:-1,0]*syx[:,0]
    G[2,1:-1,0] = p[1:-1,0]*syy[:,0]
    G[3,:,0] = 0

    return (G)
######################
## Flux Calculation ##
######################
def flux(u_in):
    
    u = zeros((4,n+2,m+1))
    u[:,1:n+1,0:m] = u_in[:,:,:]

    #Branch Cut BC
    u[:,-1,0:m] = u_in[:,0,:]
    u[:,0,0:m] = u_in[:,-1,:]
    #Airfoil BC
    #Impose this in the Flux definition
    
    #Outer BC
#    print "******************************"
#    print "** Dont forget the outer BC **"
#    print "******************************"
    
    F = f_func(u)
    G = g_func(u)

    #Initialize Fluxes
    fRight = zeros((4,n,m))
    fLeft = zeros((4,n,m))
    fUp = zeros((4,n,m))
    fDown = zeros((4,n,m))

#    print shape(fUp[:,0:n+1,1:m])
#    print shape(F[:,0:n,1:m])
#    print shape(1./2*(F[:,0:n,1:m]+F[:,0:n,2:m+1])+1./2*(G[:,0:n,1:m]+G[:,0:n,2:m+1]))

    fRight[:,0:n+1,:] = 1./2*(F[:,1:n+1,1:m+1]+F[:,2:n+2,1:m+1])*sxx[:,:]+1./2*(G[:,1:n+1,1:m+1]+G[:,2:n+2,1:m+1])*sxy[:,:]
    fLeft[:,0:n+1,:] = 1./2*(F[:,0:n,1:m+1]+F[:,1:n+1,1:m+1])*sxx[:,:]+1./2*(G[:,0:n,1:m+1]+G[:,1:n+1,1:m+1])*sxy[:,:]
    fUp[:,0:n+1,1:m] = 1./2*(F[:,0:n,1:m]+F[:,0:n,2:m+1])*sxx[:,1:m]+1./2*(G[:,0:n,1:m]+G[:,0:n,2:m+1])*sxy[:,1:m]
    fDown[:,0:n+1,1:m] = 1./2*(F[:,0:n,0:m-1]+F[:,0:n,1:m])*sxx[:,0:m-1]+1./2*(G[:,0:n,0:m-1]+G[:,0:n,1:m])*sxy[:,0:m-1]

    Flux = fRight+fLeft+fUp+fDown
    return Flux

####################
## Tau Calculator ##
####################
def tau_func(u):

    tau = zeros((n+2,m+1))
    c = zeros((n,m))
    p = zeros((n,m))

    p[:,:] = (g-1)*(u[3,0:n,0:m]-(u[1,0:n,0:m]**2+u[2,0:n,0:m]**2)/(2*u[0,0:n,0:m]))
    c[:,:] = sqrt(p[:,:]/u[0,0:n,0:m])
    tau[1:n+1,0:m] = CFL/(abs((u[1,1:n+1,0:m]/u[0,1:n+1,0:m]+c[:,:])*sxx[:,:]+(u[2,1:n+1,0:m]/u[0,1:n+1,0:m]+c[:,:])*sxy[:,:])+abs((u[1,1:n+1,0:m]/u[0,1:n+1,0:m]+c[:,:])*syx[:,:]+(u[2,1:n+1,0:m]/u[0,1:n+1,0:m]+c[:,:])*syy[:,:]))
    return(tau)

##################
## Mesh Plotter ##
##################
def meshPlotter(xs,ys,x,y):
    figure(1)
    for i in range(0,N):
        plot(xs[i,:],ys[i,:],'b')
    for j in range(0,M):
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
CFL = 1.6

print "Loading Grid Points..."
(x_mesh,y_mesh) = loadXY()

print "-> The shape of the orig girds are: [%r %s]"%(shape(x_mesh)[0],shape(x_mesh)[1])

#N,M are from shape of inputted grid
N = shape(x_mesh)[0]
M = shape(x_mesh)[1]

print "Computing Cell Centers..."
(x,y) = cellCenters(x_mesh,y_mesh)
print "-> The shape of the cell center grids are: [%r %s]"%(shape(x)[0],shape(x)[1])
print "Computing Cell Areas..."
omega = cellAreas(x_mesh,y_mesh)
print "Computing Cell Normals..."
(sxx,sxy,syx,syy) = cellNorms(x_mesh,y_mesh)

#meshPlotter(x_mesh,y_mesh,x,y)
#normalPlotter(x,y,x_mesh,y_mesh,sxx,sxy,syx,syy)
print "Initializing vectors..."

#Set IC's and BC's together assuming an initial uniform velocity field
u = zeros((4,n+2,m+1))
F = zeros((4,n+2,m+1))
G = zeros((4,n+2,m+1))
u[0,:,:] = 1.0*1000#initialize rho
u[1,:,:] = M_stream#initialize x velocity
u[2,:,:] = 0#initialize y velocity
u[3,:,:] = p_0/(g-1)+rho*(u[1,:,:]**2+u[2,:,:]**2)/2#initialize energy

a1 = 1./4
a2 = 1./3
a3 = 1./2
a4 = 1./1
tau = zeros((4,n+2,m+1))
tau[:,:] = tau_func(u)
u1 = zeros((4,n+2,m+1))
u2 = zeros((4,n+2,m+1))
u3 = zeros((4,n+2,m+1))


"""

for i in range(0,2):
    
    u1[:,:,:] = u[:,:,:] - a1*tau[:,:]*flux(u[:,:,:])
    u2[:,:,:] = u[:,:,:] - a2*tau[:,:]*flux(u1[:,:,:])
    u3[:,:,:] = u[:,:,:] - a3*tau[:,:]*flux(u2[:,:,:])
    u[:,:,:] = u[:,:,:] - a4*tau[:,:]*flux(u3[:,:,:])
    
print time.time() - start_time, "seconds"
"""
