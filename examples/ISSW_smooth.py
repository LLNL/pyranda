from __future__ import print_function
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM, pyrandaProbes


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except:
    Npts = 128

try:
    test = bool(int(sys.argv[2]))
except:
    test = False

try:
    testName = (sys.argv[3])
except:
    testName = None


## Define a mesh
L = numpy.pi * 2.0  
dim = 2
gamma = 1.4

problem = 'ISSW_smoothOperator'

Lp = L * (Npts-1.0) / Npts

from meshTest import zoomMesh_solve

dxf = 4*Lp / float(Npts) * .3
xS = zoomMesh_solve(Npts,-2.*Lp,2.*Lp,-2.,2.,1.0,dxf)

def zoomMesh(i,j,k):
    x = xS[i]
    y = xS[j]
    z = 0.0
    return x,y,z


mesh_options = {}
#mesh_options['coordsys'] = 3
#mesh_options['function'] = zoomMesh
mesh_options['periodic'] = numpy.array([False, False, False])
mesh_options['periodicGrid'] = False
mesh_options['dim'] = 3
mesh_options['x1'] = [ -2*Lp , -2*Lp  ,  0.0 ]
mesh_options['xn'] = [ 2*Lp   , 2*Lp    ,  Lp ]
mesh_options['nn'] = [ Npts, Npts ,  1  ]


# Initialize a simulation object on a mesh
ss = pyrandaSim(problem,mesh_options)
ss.addPackage( pyrandaBC(ss) )
ss.addPackage( pyrandaIBM(ss) )
ss.addPackage( pyrandaTimestep(ss) )


rho0 = 1.0
p0   = 1.0
gamma = 1.4
mach = 2.0
s0 = numpy.sqrt( p0 / rho0 * gamma )
u0 = s0 * mach
e0 = p0/(gamma-1.0) + rho0*.5*u0*u0


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  :rho:*0.0
:rho: = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
[:rx:,:ry:,:rz:] = grad( :rho: )
:v1: = :rx:*:gx: + :ry:*:gy:
"""
eom = eom.replace('u0',str(u0)).replace('p0',str(p0)).replace('rho0',str(rho0))


# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = """
rad = sqrt( meshx**2  +  meshy**2 ) 
theta = numpy.arctan2(meshy,meshx)
:phi: = rad - pi
#:rho: = 1.0 + 3d()
#:rho: = meshx
:rho: = sin( theta )
[:gx:,:gy:,:gz:] = grad( :phi: )
"""

# Set the initial conditions
ss.setIC(ic)
    
# Length scale for art. viscosity
# Initialize variables
x = ss.mesh.coords[0].data
y = ss.mesh.coords[1].data
z = ss.mesh.coords[2].data

# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
tt = .2 #

# Start time loop
cnt = 1
viz_freq = 10
pvar = 'umag'


# Make a probe set for diagnostics
# Make a probe set for diagnostics
nProbes = Npts
theta = numpy.linspace(0,2.0*numpy.pi,nProbes+1)[:-1]
R0 = numpy.pi
x = R0 * numpy.cos( theta )
y = R0 * numpy.sin( theta )
probes = pyrandaProbes(ss,x=x,y=y,z=None)



while tt > time:
    
    # Update the EOM and get next dt
    dt = 0.01
    time = ss.rk4(time,dt)
    
    # Print some output
    #ss.iprint("%s -- %s --- %f" % (cnt,time,dt)  )

    prb = probes.get('v1')
    print("Var: %s    Max: %s " % (prb.var(),prb.max()   ))
    
    
    if viz and (not test):
        if (cnt%viz_freq == 1) :
            ss.write()
        
            #ss.plot.figure(1)
            #ss.plot.clf()            
            #ss.plot.contourf('v1' ,64 , cmap=cm.jet)
            #ss.plot.contour('phi',[0.0])

            #ss.plot.figure(2)
            
    cnt += 1


ss.plot.figure(1)
ss.plot.clf()            
ss.plot.contourf('v1' ,64 , cmap=cm.jet)
ss.plot.contour('phi',[0.0])


ss.plot.figure(2)
prb = probes.get('v1')
plt.plot( prb )
plt.show()
