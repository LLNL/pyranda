from __future__ import print_function
import sys,os
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM

# Try to get args
try:
    Npts = int(sys.argv[1])
except:
    Npts = 27

try:
    test = bool(sys.argv[2])
except:
    test = False

## Define a mesh
L = numpy.pi * 2.0  
Lp = L * (Npts-1.0) / Npts


mesh_options = {}
mesh_options['coordsys'] = 0
mesh_options['periodic'] = numpy.array([False, True, True])
mesh_options['dim'] = 1
mesh_options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
mesh_options['xn'] = [ Lp   , Lp    ,  Lp ]
mesh_options['nn'] = [ Npts, 1 ,  1  ]

# Initialize a simulation object on a mesh
ss = pyrandaSim('heat_equation',mesh_options)

ss.addPackage( pyrandaBC(ss) )
ss.addPackage( pyrandaIBM(ss) )
ss.addPackage( pyrandaTimestep(ss) )


# Define the equations of motion
eom = """
# Heat equation
ddt(:T:)  = lap(:T:)
# Boundary condition
:T: = ibmC( :T: , :phi:, [:gx:,:gy:,:gz:] , :Tsol2:)
"""
ss.EOM(eom)

# Initialize variables
x1 = numpy.pi/4.0
x2 = numpy.pi*7.0/4.0
T1 = 1.0
T2 = 2.0

ic = """
x1 = pi/4.0
x2 = pi*7.0/4.0
T1 = 1.0
T2 = 2.0
:Tsol: = T1 + (T2-T1)*(x1-meshx)/(x1-x2)
:Tsol2: = where(meshx<pi,T1,T2)
#:T: = where(meshx<x1,T1,:Tsol: )
#:T: = where(meshx>x2,T2,:T:    )
:T: = :Tsol:*1.0
:phi:  = numpy.minimum( meshx - x1 , x2 - meshx )
[:gx:,:gy:,:gz:] = grad( :phi: )
"""
ss.setIC(ic)

x  = ss.mesh.coords[0].data
xx =  ss.PyMPI.zbar( x )

# Time step size
dt_max = L / ss.mesh.nn[0] * .005
tt = dt_max * 2000 

rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 20
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

# Main time loop for physics
dt = dt_max
cnt = 1
time = 0.0
viz = True
while tt > time:

    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )

    # Plot animation of advection
    anl = ss.PyMPI.zbar( ss.variables['Tsol'].data )
    v = ss.PyMPI.zbar( ss.variables['T'].data )
    
    error = numpy.abs( anl - v )
    E2 = numpy.mean(numpy.abs(error[:,0]))

    if not test:
        ss.iprint("%s -- %s -- %s" % (cnt,time,E2)  )
        
    cnt += 1
    if viz:
        if (ss.PyMPI.master and (cnt%50 == 0)) and (not test):
            plt.figure(1)
            plt.clf()

            LW = 2.5
            plt.plot(xx[:,0],v[:,0],'bo-',linewidth=LW)
            #plt.plot(xx[:,0],anl[:,0],'k.-',linewidth=LW)
            plt.plot( [x1,x1], [0.5,2.5] , 'k--' )
            plt.plot( [x2,x2], [0.5,2.5] , 'k--' )
            plt.plot( [0,x1],[T1,T1], 'r-',linewidth=LW)
            plt.plot( [x2,2*numpy.pi],[T2,T2], 'r-',linewidth=LW)
            plt.pause(.001)



# Make pretty
floc = "/Users/olson45/Documents/Conference Travel/ISSW32_2019/ISSW32_bjo/figures"

plt.xlabel("$x$") 
plt.ylabel("$T(x)$")
plt.tight_layout(pad=1.0)
#plt.savefig(os.path.join(floc,'heat_1D.pdf'))


foo = """
17  0.000376208943829
27  1.50894351884e-05
54  1.68226618088e-08
108 2.02899425656e-13
217 1.9937968327e-15

res = [ [17  ,0.000376208943829],
        [27  ,1.50894351884e-05],
        [54  ,1.68226618088e-08],
        [108 ,2.02899425656e-13],
        [217 ,1.9937968327e-15]]

import numpy as npy
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os,sys

res = npy.array( res )

rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 20
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
floc = "/Users/olson45/Documents/Conference Travel/ISSW32_2019/ISSW32_bjo/figures"


LW = 2.5

rate = 10
invdx = res[:,0]/(2.0*npy.pi)
scl = invdx[0]**(-rate)/res[0,1]
plt.loglog( invdx , res[:,1], 'bo', linewidth=LW)
plt.loglog( invdx , 8.0/scl * invdx**(-rate), 'k--')


plt.xlabel(r"$\frac{1}{\Delta x}$") 
plt.ylabel(r"$L_2$ Error")
plt.tight_layout(pad=1.0)
plt.savefig(os.path.join(floc,'heat_1D_conv.pdf'))

plt.show()



"""
