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
    airfoilData[:,1] = airfoilData[:,1]-0.5

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
#    for i in range(0,n):
#        x[i,0] = airfoil[i,1]
#        y[i,0] = airfoil[i,2]
    x[:,0] = airfoil[:,1]
    y[:,0] = airfoil[:,2]

    #Set outer BC to circle with a radius of 10 cord lenths
    dTheta = 2*pi/(n-1) #difference in angle between points
    for i in range(0,n-1):
        x[i,m-1] = r*cos(dTheta*i)
        y[i,m-1] = -r*sin(dTheta*i)
 
    #Set location of "cut"
#    dx = (r - airfoil[0,1])/m #spacing from tail to boundary

    tip = airfoil[0,1]
    dx = float(r-tip)/(m-1)
    for j in range(1,m):
        x[0,j] = tip+dx*j
        y[0,j] = 0
        x[n-1,j] = tip+dx*j
        y[n-1,j] = 0

    return (x,y)

 
################
## BC Plotter ##
################
def bcPlotter(x,y):
    figure(2)
    title('Mapping of Boundaries to Physical Domain')
    xlabel('x')
    ylabel('y')
    plot(x[:,0],y[:,0],'r')
    plot(x[:,m-1],y[:,m-1],'g')
    plot(x[0,:],y[0,:],'b')
    plot(x[n-1,:],y[n-1,:],'k--')
    legend(('Bottom Boundary','Top Boundary','Left Boundary','Right Boundary'))
    
    return

##########################
## Initial Mesh Plotter ##
##########################
def initMeshPlotter(x,y):
    figure(3)
    title('Mapping of Initial Mesh')
    xlabel('x')
    ylabel('y')
    plot(x[:,:],y[:,:],'b.')
    plot(x[:,0],y[:,0],'r.')
    plot(x[0,0],y[0,0],'g.')
    axis([-10.5,10.5,-10.5,10.5]) 
#    legend(('Transfinite Interpolation','NACA 0012 Airfoil'))
    savefig('initMesh')
#    axis([-0.5,1.5,-0.5,0.5]) 
#    savefig('initMeshZoom')
    return

##################
## Mesh Plotter ##
##################
def meshPlotter(x,y):
    figure(4)
    title('Mapping of Computaitonal Domain to Physical Domain')
    xlabel('x')
    ylabel('y')
    plot(x[:,:],y[:,:],'b.')
    plot(x[:,0],y[:,0],'r')
    savefig('meshPlot')
    return

######################
## Residual Plotter ##
######################
def resPlotter(res):
    figure(5)
    title('residual')
    xlabel('n')
    ylabel('Residual')
    plot(res)
    savefig('resPlot')
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

#            L1j = float(j-m)/(0-m)
#            L2j = float(j-0)/(m-0)
#            L1i = float(i-n)/(0-n)
#            L2i = float(i-0)/(n-0)


            L1j = float(j-(m-1))/(0-(m-1))
            L2j = float(j-0)/((m-1)-0)
            L1i = float(i-(n-1))/(0-(n-1))
            L2i = float(i-0)/((n-1)-0)


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

################################
## Successive Over Relaxation ##
################################
def SOR(Rmin,omega,x,y):
    res = [] #array of residuals
    R = Rmin*2
    while R>=Rmin:
        R = 0
       
        for j in range(1,m-1):
            for i in range(1,n-1):
             
                
                xold[i,j] = x[i,j]
                yold[i,j] = y[i,j]
   
                alpha = 0.25*((x[i,j+1]-x[i,j-1])**2+(y[i,j+1]-y[i,j-1])**2)
                beta =  0.25*((x[i+1,j]-x[i-1,j])*(x[i,j+1]-x[i,j-1])+(y[i+1,j]-y[i-1,j])*(y[i,j+1]-y[i,j-1]))
                gamma = 0.25*((x[i+1,j]-x[i-1,j])**2+(y[i+1,j]-y[i-1,j])**2)
                delta = 1./16*((x[i+1,j]-x[i-1,j])*(y[i,j+1]-y[i,j-1])-(x[i,j+1]-x[i,j-1])*(y[i+1,j]-y[i-1,j]))**2

                P = 0
                Q = 0
                #SOR step
                x[i,j] = 1./(2*(alpha+gamma))*(alpha*(x[i-1,j]+x[i+1,j])-0.5*beta*(x[i+1,j+1]-x[i-1,j+1]-x[i+1,j-1]+x[i-1,j-1])+gamma*(x[i,j-1]+x[i,j+1])+0.5*delta*P*(x[i+1,j]-x[i-1,j])+0.5*delta*Q*(x[i,j+1]-x[i,j-1]))
                y[i,j] = 1./(2*(alpha+gamma))*(alpha*(y[i-1,j]+y[i+1,j])-0.5*beta*(y[i+1,j+1]-y[i-1,j+1]-y[i+1,j-1]+y[i-1,j-1])+gamma*(y[i,j-1]+y[i,j+1])+0.5*delta*P*(y[i+1,j]-y[i-1,j])+0.5*delta*Q*(y[i,j+1]-y[i,j-1]))
                                
                #Relaxation step
                x[i,j] = omega*x[i,j]+(1-omega)*xold[i,j]
                y[i,j] = omega*y[i,j]+(1-omega)*yold[i,j]
                



        for j in range(1,m-1):
            for i in range(1,n-1):
                
                alpha = 0.25*((x[i,j+1]-x[i,j-1])**2+(y[i,j+1]-y[i,j-1])**2)
                beta =  0.25*((x[i+1,j]-x[i-1,j])*(x[i,j+1]-x[i,j-1])+(y[i+1,j]-y[i-1,j])*(y[i,j+1]-y[i,j-1]))
                gamma = 0.25*((x[i+1,j]-x[i-1,j])**2+(y[i+1,j]-y[i-1,j])**2)
                delta = 1./16*((x[i+1,j]-x[i-1,j])*(y[i,j+1]-y[i,j-1])-(x[i,j+1]-x[i,j-1])*(y[i+1,j]-y[i-1,j]))**2

                #Compute residual
                Rx = (alpha*(x[i-1,j]-2*xold[i,j]+x[i+1,j])-0.5*beta*(x[i+1,j+1]-x[i-1,j+1]-x[i+1,j-1]+x[i-1,j-1])+gamma*(x[i,j-1]-2*xold[i,j]+x[i,j+1])+0.5*delta*P*(x[i+1,j]-x[i-1,j])+0.5*delta*Q*(x[i,j+1]-x[i,j-1]))
                Ry = (alpha*(y[i-1,j]-2*yold[i,j]+y[i+1,j])-0.5*beta*(y[i+1,j+1]-y[i-1,j+1]-y[i+1,j-1]+y[i-1,j-1])+gamma*(y[i,j-1]-2*yold[i,j]+y[i,j+1])+0.5*delta*P*(y[i+1,j]-y[i-1,j])+0.5*delta*Q*(y[i,j+1]-y[i,j-1]))

#                Rx = x[i,j] - 2*(alpha+gamma)*x[i,j]
#                Ry = y[i,j] - 2*(alpha+gamma)*y[i,j]


                R = (abs(Rx)+abs(Ry))/2+R
            
        R = R/((n-1)*(m-1))
        res.append(R)

        print len(res),R
        if len(res)>1300:
            break

#    print Rx, Ry
    return (x,y,res)




##################
## Main Program ##
##################

close ('all')
n = 129
m = 65
r = 10 #(10 cord lengths) x (cord lenght = 1) 
Rmin = 10**(-6)
omega = 1.5
airfoil = loadAirfoil()
x,y = initBC()
#bcPlotter(x,y)
#plotAirfoil(airfoil)

x,y = initMesh(x,y)
initMeshPlotter(x,y)
xold = x
yold = y
x,y,res = SOR(Rmin,omega,x,y)
meshPlotter(x,y)
resPlotter(res)
show()
