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
    sxx[0:n,0:m] = (ys[1:N,1:M]-ys[1:N,0:M-1])
    sxy[0:n,0:m] = -(xs[1:N,1:M]-xs[1:N,0:M-1])
    syx[0:n,0:m] = -sxy[:,:]
    syy[0:n,0:m] = sxx[:,:]

    return (sxx,sxy,syx,syy)

################
## F Function ##
################
def f_func(u):
    N = shape(u)[1]
    M = shape(u)[2]
    F = zeros((4,N,M))
    p = zeros((N,M))
    
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
def flux(u):
    
    #Only need to compute F and G once
    F = f_func(u)
    G = g_func(u)

    #Initialize Fluxes
    fRight = zeros((4,n,m))
    fLeft = zeros((4,n,m))
    fUp = zeros((4,n,m))
    fDown = zeros((4,n,m))

    #Initialize Dissapative Terms
    dRight = zeros((4,n,m))
    dLeft = zeros((4,n,m))
    dUp = zeros((4,n,m))
    dDown = zeros((4,n,m))
    maxNu = zeros((n,m))

    #Compute flux for each direction individually without dissapation
    fRight[:,:,:] = 1./2*((F[:,1:n+1,1:m+1]+F[:,2:n+2,1:m+1])*sxx[:,:]+(G[:,1:n+1,1:m+1]+G[:,2:n+2,1:m+1])*sxy[:,:])
    fLeft[:,:,:] =  1./2*((F[:,0:n,1:m+1]+F[:,1:n+1,1:m+1])*sxx[:,:]+(G[:,0:n,1:m+1]+G[:,1:n+1,1:m+1])*sxy[:,:])
    fUp[:,:,:] =    1./2*((F[:,1:n+1,1:m+1]+F[:,1:n+1,2:m+2])*sxx[:,:]+(G[:,1:n+1,1:m+1]+G[:,1:n+1,2:m+2])*sxy[:,:])
    fDown[:,:,:] =  1./2*((F[:,1:n+1,1:m+1]+F[:,1:n+1,0:m])*sxx[:,:]+(G[:,1:n+1,1:m+1]+G[:,1:n+1,0:m])*sxy[:,:])

    #Compute dissapation for each direction

    (eps2i,eps2j,eps4i,eps4j) = eps_calc(u)
    dRight[:,:,:] = eps2i[:,:]*(u[:,2:n+2,1:m+1]-u[:,1:n+1,1:m+1]) - eps4i[:,:]*(u[:,3:mod(n+3,n)+1,1:m+1]-3*u[:,2:n+2,1:m+1]+3*u[:,1:n+1,1:m+1]-u[:,0:n,1:m+1])
    dLeft[:,1:n,1:m] = -dRight[:,0:n-1,1:m]
    dUp[:,:,:] = eps2j[:,:]*(u[:,1:n+1,2:m+2]-u[:,1:n+1,1:m+1]) - eps4i[:,:]*(u[:,1:n+1,3:mod(m+3,m)+1]-3*u[:,1:n+1,2:m+2]+3*u[:,1:n+1,1:m+1]-u[:,1:n+1,0:m])
    dDown[:,1:n,1:m] = -dUp[:,1:n,0:m-1]

    Flux = fRight+fLeft+fUp+fDown-(dRight+dLeft+dUp+dDown)
    return Flux

###################
## Nu Calculator ##
##################
def nu_max_func(u):
    #This function computes the maximum nu required for the dissapation calculation and return an n x m vector of values
    
    nui = zeros((n,m))
    nuj = zeros((n,m))
    maxNui = zeros((n,m))
    maxNuj = zeros((n,m))
    p = zeros((n+2,m+2))

    p[:,:] = (g-1)*(u[3,:,:]-(u[1,:,:]**2+u[2,:,:]**2)/(2*u[0,:,:]))    

    nui[:,:] = abs((p[0:n,1:m+1]-2*p[1:n+1,1:m+1]+p[2:n+2,1:m+1])/(p[0:n,1:m+1]+2*p[1:n+1,1:m+1]+p[2:n+2,1:m+1]))
    nuj[:,:] = abs((p[1:n+1,0:m]-2*p[1:n+1,1:m+1]+p[1:n+1,2:m+2])/(p[1:n+1,0:m]+2*p[1:n+1,1:m+1]+p[1:n+1,2:m+2]))

    maxNui[1:n-2,:] = maximum(maximum(nui[0:n-3,:],nui[1:n-2,:]),maximum(nui[2:n-1,:],nui[3:n,:]))
    maxNui[n-2,:] = maximum(nui[n-3,:],maximum(nui[n-2,:],nui[n-1,:]))
    maxNui[n-1,:] = maximum(nui[n-2,:],nui[n-1,:])

    maxNuj[:,1:m-2] = maximum(maximum(nuj[:,0:m-3],nuj[:,1:m-2]),maximum(nuj[:,2:m-1],nuj[:,3:m]))
    maxNuj[:,m-2] = maximum(nuj[:,m-3],maximum(nuj[:,m-2],nuj[:,m-1]))
    maxNuj[:,m-1] = maximum(nuj[:,m-2],nui[:,m-1])
    return (maxNui,maxNuj)

########################
## Epsilon Calculator ##
########################
def eps_calc(u):
    eps2i = zeros((n,m))
    eps4i = zeros((n,m))
    eps2j = zeros((n,m))
    eps4j = zeros((n,m))
    p = zeros((n,m))
    c = zeros((n,m))

    k2 = 1./4
    k4 = 1./256

    u_vel = 1./2*(u[1,1:n+1,1:m+1]/u[0,1:n+1,1:m+1]+u[1,2:n+2,1:m+1]/u[0,2:n+2,1:m+1])
    v_vel = 1./2*(u[2,1:n+1,1:m+1]/u[0,1:n+1,1:m+1]+u[2,1:n+1,2:m+2]/u[0,1:n+1,2:m+2])

    p[:,:] = (g-1)*(u[3,1:n+1,1:m+1]-(u_vel[:,:]**2+v_vel[:,:]**2)/(2*u[0,1:n+1,1:m+1]))
    print "The max and min of p/p_0 are", p.max()/p_0, p.min()/p_0
                
    c[:,:] = sqrt(p[:,:]/u[0,1:n+1,1:m+1])
    
    print "The max and min  of c are", c.max(), c.min()
    #Compute Max Nu values to be used in epsilon calculations
    (maxNui,maxNuj) = nu_max_func(u)

    eps2i[:,:] = 1./2*k2*(u_vel[:,:]*sxx[:,:]+v_vel*sxy[:,:]+c[:,:]*sqrt(sxx[:,:]**2+sxy[:,:]**2))*maxNui[:,:]
    eps2j[:,:] = 1./2*k2*(u_vel[:,:]*syx[:,:]+v_vel*syy[:,:]+c[:,:]*sqrt(syx[:,:]**2+syy[:,:]**2))*maxNuj[:,:]
    
    eps4i[:,:] = maximum(zeros((n,m)),1./2*k4*(u_vel[:,:]*sxx[:,:]+v_vel*sxy[:,:]+c[:,:]*sqrt(sxx[:,:]**2+sxy[:,:]**2))-eps2i[:,:])
    eps4j[:,:] = maximum(zeros((n,m)),1./2*k4*(u_vel[:,:]*syx[:,:]+v_vel*syy[:,:]+c[:,:]*sqrt(syx[:,:]**2+syy[:,:]**2))-eps2j[:,:])

    return(eps2i,eps2j,eps4i,eps4j)

####################
## Tau Calculator ##
####################
def tau_func(u):
    
    #This function computes the time step in each cell normalized by cell area (not needed explicitly) and puts zeros in the ghost cell locations
    tau = zeros((n,m))
    c = zeros((n,m))
    p = zeros((n,m))
    dsxi = zeros((n,m))
    dsxj = zeros((n,m))
    dsyi = zeros((n,m))
    dsyj = zeros((n,m))
       
    dsxi[:,:] = 1./2*(sxx[0:n,:]+sxx[range(1,n+1).append(0),:])
    dsxj[:,1:m] = 1./2*(sxy[:,0:m-1]+sxy[:,1:m])
    dsxj[:,0] = 0

    dsyi[:,:] = 1./2*(syx[0:n,:]+syx[range(1,n+1).append(0),:])
    dsyj[:,1:m] = 1./2*(syy[:,0:m-1]+syy[:,1:m])
    dsyj[:,0] = 0

    p[:,:] = (g-1)*(u[3,1:n+1,1:m+1]-(u[1,1:n+1,1:m+1]**2+u[2,1:n+1,1:m+1]**2)/(2*u[0,1:n+1,1:m+1]))
    c[:,:] = sqrt(p[:,:]/u[0,1:n+1,1:m+1])

    print shape(dsxi[:,:])
    a = abs((abs(u[1,1:n+1,1:m+1]/u[0,1:n+1,1:m+1])+c[:,:])*dsxi[:,:])
    b = abs((abs(u[2,1:n+1,1:m+1]/u[0,1:n+1,1:m+1])+c[:,:])*dsxj[:,:])
    c = abs((abs(u[1,1:n+1,1:m+1]/u[0,1:n+1,1:m+1])+c[:,:])*dsyi[:,:])
    d = abs((abs(u[2,1:n+1,1:m+1]/u[0,1:n+1,1:m+1])+c[:,:])*dsyj[:,:])
    tau[:,:] = CFL/(a+b+c+d)

    return(tau)

#################################
## Boundary Condition Enforcer ##
#################################

def bc_enforce(u):
    #Note that Airfoil BC is enforced inside of flux functions and hence are not addressed by this function

    #Branch Cut BC
    u[:,-1,:] = u[:,1,:]
    u[:,0,:] = u[:,-2,:]
    
    #Outer BC
    pInt = zeros(n)
    pInt[:] = (g-1)*(u[3,1:n+1,-2]-(u[1,1:n+1,-2]**2+u[2,1:n+1,-2]**2)/(2*u[0,1:n+1,-2]))
    for i in range(1,n):
        #Inflow Condition
        if syx[i,-1] < 0:
            uni = u[1,i,-2]/u[0,i,-2]*syx[i,-2]
            rnb = M_stream*syx[i,-2]+2*sqrt(p_0/rho)/(g-1)
            rni = uni-2*sqrt(pInt[i]/u[0,i,-2])/(g-1)
            cb =  1./2*(rnb-rni)*(g-1)/4
            unb = 1./2*(rnb+rni)
            utb = M_stream*syx[i,-2]

            #Rho Boundary
            u[0,i,-1] = (cb**2/p_0*rho**g)**(1/(g-1))
            #pressure on Boundary
            pb = cb**2*u[0,i,-2]
            #rho*u Boundary
            u[1,i,-1] = u[0,i,-2]*(unb*syx[i,-2]+utb*syx[i,-2])
            #rho*v Boundary
            u[2,i,-1] = u[0,i,-2]*(unb*syy[i,-2]+utb*syy[i,-2])
            #rho*E Boundary
            u[3,i,-1] = pInt[i]/(g-1)+u[0,i,-2]*((u[1,i,-2]/u[0,i,-2])**2+(u[2,i,-2]/u[0,i,-2])**2)/2

        #Outflow Condition
        else:
            uni = u[1,i,-1]/u[0,i,-1]*syx[i,-1]
            rnb =-u[1,i,-1]/u[0,i,-1]*syx[i,-1]+2*sqrt(pInt[i]/u[0,i,-1])/(g-1)
            rni = uni-2*sqrt(pInt[i]/u[0,i,-1])/(g-1)
            cb =  1./2*(rnb-rni)*(4-1)/4
            unb = 1./2*(rnb+rni)
            utb = M_stream*syx[i,-1]

            #Rho Boundary
            u[0,i,-1] = (cb**2/pInt[i]*u[0,i,-1]**g)**(1/(g-1))
            #pressure on Boundary
            pb = cb**2*u[0,i,-1]
            #rho*u Boundary
            u[1,i,-1] = u[0,i,-1]*(unb*syx[i,-1]+utb*syx[i,-1])
             #rho*v Boundary
            u[2,i,-1] = u[0,i,-1]*(unb*syy[i,-1]+utb*syy[i,-1])
            #rho*E Boundary
            u[3,i,-1] = pInt[i]/(g-1)+u[0,i,-1]*((u[1,i,-1]/u[0,i,-1])**2+(u[2,i,-1]/u[0,i,-1])**2)/2
   
    return u



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
    plot(x,y,'r.')
    quiver(x,y,sxx*scale,sxy*scale)
    quiver(x,y,syx*scale,syy*scale)
    axis([-10.1,-8.5,-1.3,1.3])
    title('Cell Centers')
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
    plot(x,y,'r.')
    quiver(x,y,x_vec,y_vec)
    title('Velocity Directions')
    xlabel('x')
    ylabel('y')

    return


####################
## Plotting Tools ##
####################
def plotResults(u):
    figure(3)
    title('Velocity') 
    contourf(x,y,sqrt((u[1,1:n+1,1:m+1]/u[0,1:n+1,1:m+1])**2+(u[2,1:n+1,1:m+1]/u1[0,1:n+1,1:m+1])**2))
    colorbar()
#    savefig('velocity_u1.png')
    show()
    return

###################
## Initial Plots ##
###################
def initPlots(u,time):
    figure(1)
    contourf(x,y,u[0,1:n+1,1:m+1])
    colorbar()
    title('%s Density'%time)
    xlabel('x')
    ylabel('y')
    savefig('%s_rho.png'%time)
    axis([-0.75, 0.75, -0.75, 0.75])
    savefig('%s_rho_zoom.png'%time)

    figure(2)
    velocity = sqrt((u[1,1:n+1,1:m+1]**2+u[1,1:n+1,1:m+1]**2)/u[0,1:n+1,1:m+1])
    contourf(x,y,velocity)
    colorbar()
    title('%s Velocity Magnitude'%time)
    xlabel('x')
    ylabel('y')
    savefig('%s_vel.png'%time)
    axis([-0.75, 0.75, -0.75, 0.75])
    savefig('%s_vel_zoom.png'%time)

    figure(3)
    p = (g-1)*(u[3,1:n+1,1:m+1]-(u[1,1:n+1,1:m+1]**2+u[2,1:n+1,1:m+1]**2)/(2*u[0,1:n+1,1:m+1]))
    contourf(x,y,p)
    colorbar()
    title('%s Pressure'%time)
    xlabel('x')
    ylabel('y')
    savefig('%s_pressure.png'%time)
    axis([-0.75, 0.75, -0.75, 0.75])
    savefig('%s_pressure_zoom.png'%time)

    figure(4)
    c = sqrt(p/u[0,1:n+1,1:m+1])
    contourf(x,y,velocity/c)
    colorbar()
    title('%s Mach Number'%time)
    xlabel('x')
    ylabel('y')
    savefig('%s_mach.png'%time)
    axis([-0.75, 0.75, -0.75, 0.75])
    savefig('%s_mach_zoom.png'%time)

    figure(5)
    contourf(x,y,p/u[0,1:n+1,1:m+1]**g)
    colorbar()
    title('%s Entropy'%time)
    xlabel('x')
    ylabel('y')
    savefig('%s_entropy.png'%time)
    axis([-0.75, 0.75, -0.75, 0.75])
    savefig('%s_entropy_zoom.png'%time)

    figure(6)
    contourf(x,y,(u[3,1:n+1,1:m+1]+p)/(p_0*u[0,1:n+1,1:m+1]))
    colorbar()
    title('%s Enthalpy'%time)
    xlabel('x')
    ylabel('y')
    savefig('%s_enthalpy.png'%time)
    axis([-0.75, 0.75, -0.75, 0.75])
    savefig('%s_enthalpy_zoom.png'%time)

    figure(7)
    contourf(x,y,(p-p_0)/(1./2*rho*M_stream**2))
    colorbar()
    title('%s Pressure Coefficient'%time)
    xlabel('x')
    ylabel('y')
    savefig('%s_press_coeff.png'%time)
    axis([-0.75, 0.75, -0.75, 0.75])
    savefig('%s_press_coeff_zoom.png'%time)

    figure(8)
    contourf(x,y,p/p_0)
    colorbar()
    title('%s Normalized Pressure'%time)
    xlabel('x')
    ylabel('y')
    savefig('%s_norm_press.png'%time)
    axis([-0.75, 0.75, -0.75, 0.75])
    savefig('%s_norm_press_zoom.png'%time)
    
       
    return
##################
## Main Program ##
##################
seterr('raise')
close('all')
start_time = time.time()
n = 128
m = 64
g = 1.4
p_0 = 1.*10**5
rho = 1.
M_stream = 0.85*sqrt(10**5/rho)
CFL = 0.2

print "Loading Grid Points..."
(x_mesh,y_mesh) = loadXY()
#N,M are from shape of inputted grid
N = shape(x_mesh)[0]
M = shape(x_mesh)[1]
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
u = zeros((4,n+2,m+2))
tau = zeros((4,n,m))
u[0,:,:] = 1.0#initialize rho
u[1,:,:] = M_stream#initialize x velocity
u[2,:,:] = 0#initialize y velocity
u[3,:,:] = p_0/(g-1)+rho*(M_stream**2)/2#initialize energy

#Ruga-Kutta Coeffs
a1 = 1./4
a2 = 1./3
a3 = 1./2
a4 = 1./1


u1 = zeros((4,n+2,m+2))
u2 = zeros((4,n+2,m+2))
u3 = zeros((4,n+2,m+2))

initPlots(u,'Initial')

for i in range(0,1):
    print i
    tau[:,:] = tau_func(u)
    print "The max and min of tau are", tau.max(), tau.min()
    
#    u = bc_enforce(u)
    u1[:,:,:] = u[:,:,:]
    u2[:,:,:] = u[:,:,:]
    u3[:,:,:] = u[:,:,:]
    u1[:,1:n+1,1:m+1] = u[:,1:n+1,1:m+1] - a1*tau[:,:]*flux(u[:,:,:])
    
    u1 = bc_enforce(u1)

    u2[:,1:n+1,1:m+1] = u[:,1:n+1,1:m+1] - a2*tau[:,:]*flux(u1[:,:,:])
    u2 = bc_enforce(u2)

    u3[:,1:n+1,1:m+1] = u[:,1:n+1,1:m+1] - a3*tau[:,:]*flux(u2[:,:,:])
    u3 = bc_enforce(u3)

    u[:,1:n+1,1:m+1] = u[:,1:n+1,1:m+1] - a4*tau[:,:]*flux(u3[:,:,:])

    print "The max and min  of u are", (u[1,:,:]/u[0,:,:]).max(), (u[1,:,:]/u[0,:,:]).min()
    print "The max and min  of rhoE are", (u[0,:,:]).max(), (u[3,:,:]).min()
   
print time.time() - start_time, "seconds"

