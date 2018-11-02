from __future__ import print_function
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except:
    Npts = 64

try:
    test = bool(int(sys.argv[2]))
except:
    test = False

try:
    testName = (sys.argv[3])
except:
    testName = None


## Define a mesh
#Npts = 32
L = numpy.pi * 2.0  
dim = 2
gamma = 1.4

problem = 'cylinder_curvilinear2'

Lp = L * (Npts-1.0) / Npts

Ri = 1.0
Rf = 5.0

NX = Npts
NY = Npts

def cylMesh(i,j,k):
    theta =  float(i) / float(NX) * 2.0 * numpy.pi
    r = float(j) / float(NY-1) * (Rf-Ri) + Ri
    x = r * numpy.cos( theta )
    y = r * numpy.sin( theta )
    z = 0.0
    return x,y,z


mesh_options = {}
mesh_options['coordsys'] = 3
mesh_options['function'] = cylMesh
mesh_options['periodic'] = numpy.array([True, False, True])
mesh_options['dim'] = 3
mesh_options['x1'] = [ -2*Lp , -2*Lp  ,  0.0 ]
mesh_options['xn'] = [ 2*Lp   , 2*Lp    ,  Lp ]
mesh_options['nn'] = [ NX, NY ,  1  ]


# Initialize a simulation object on a mesh
ss = pyrandaSim(problem,mesh_options)
ss.addPackage( pyrandaBC(ss) )
ss.addPackage( pyrandaIBM(ss) )
ss.addPackage( pyrandaTimestep(ss) )


rho0 = 1.0
p0   = 1.0
gamma = 1.4
mach = 1.2
s0 = numpy.sqrt( p0 / rho0 * gamma )
u0 = s0 * mach
e0 = p0/(gamma-1.0) + rho0*.5*u0*u0


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -div(:rho:*:u:,  :rho:*:v:)
ddt(:rhou:) =  -div(:rhou:*:u: + :p: - :tau:, :rhou:*:v:)
ddt(:rhov:) =  -div(:rhov:*:u:, :rhov:*:v: + :p: - :tau:)
ddt(:Et:)   =  -div( (:Et: + :p: - :tau:)*:u: - :tx:*:kappa:, (:Et: + :p: - :tau:)*:v: - :ty:*:kappa: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:,:v:)
:beta:      =  gbar( ring(:div:) * :rho: ) * 7.0e-2
:tau:       = :beta:*:div:
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
bc.extrap(['u','v','rho','p'],['yn'])
bc.extrap(['rho','p'],['y1'])
bc.slip([ ['u','v']  ],['y1'])
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 2.2* dt.diff(:beta:,:rho:)
:dt: = numpy.minimum(:dt:,:dtB:)
:umag: = sqrt( :u:*:u: + :v:*:v: )
"""
eom = eom.replace('u0',str(u0)).replace('p0',str(p0)).replace('rho0',str(rho0))


# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = """
:gamma: = 1.4
:R: = 1.0
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
#rad = sqrt( (meshx-numpy.pi)**2  +  (meshy-numpy.pi)**2 ) 
rad = sqrt( meshx**2  +  meshy**2 ) 
:phi: = rad - numpy.pi/4.0
:rho: = 1.0 + 3d()
:p:  =  1.0 + 3d() #exp( -(meshx-1.5)**2/.25**2)*.1
#:u: = where( :phi:>0.5, mach * sqrt( :p: / :rho: * :gamma:) , 0.0 )
:u: = mach * sqrt( :p: / :rho: * :gamma:)
#:u: = gbar( gbar( :u: ) )
:v: = 0.0 + 3d()
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
[:gx:,:gy:,:gz:] = grad( :phi: )
:gx: = gbar( :gx: )
:gy: = gbar( :gy: )
"""
ic = ic.replace('mach',str(mach))

# Set the initial conditions
ss.setIC(ic)
    
# Length scale for art. viscosity
# Initialize variables
x = ss.mesh.coords[0]
y = ss.mesh.coords[1]
z = ss.mesh.coords[2]
dx = (x[1,0,0] - x[0,0,0])
#ss.variables['dx2'].data += dx**2


# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
tt = 3.0 #

# Mesh for viz on master
xx   =  ss.PyMPI.zbar( x )
yy   =  ss.PyMPI.zbar( y )

# Start time loop
cnt = 1
viz_freq = 10
pvar = 'umag'

#import pdb
#pdb.set_trace()
CFL = 0.2
dt = ss.variables['dt'].data * CFL




v = ss.PyMPI.zbar( ss.variables[pvar].data )
phi = ss.PyMPI.zbar( ss.variables['phi'].data )
if (ss.PyMPI.master and (not test) ):


    plt.figure(2)
    plt.clf()            
    plt.contourf( xx,yy,v ,64 , cmap=cm.jet)
    plt.contour( xx,yy,phi,[0.0])
    plt.plot(xx, yy, 'k-', lw=0.5, alpha=0.5)
    plt.plot(xx.T, yy.T, 'k-', lw=0.5, alpha=0.5)
    plt.title(pvar)
    plt.pause(.001)



while tt > time:
    
    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = ss.variables['dt'].data * CFL
    dt = min(dt, (tt - time) )

    
    # Print some output
    ss.iprint("%s -- %s --- %f" % (cnt,time,dt)  )
    cnt += 1
    if viz and (not test):
        v = ss.PyMPI.zbar( ss.variables[pvar].data )
        phi = ss.PyMPI.zbar( ss.variables['phi'].data )
        if (cnt%viz_freq == 1) :
            ss.write(['beta','u','v','rho','umag'])
        
        if (ss.PyMPI.master) and (cnt%viz_freq == 1) :#or True:

            
            
            plt.figure(2)
            plt.clf()            
            plt.contourf( xx,yy,v ,64 , cmap=cm.jet)
            plt.contour( xx,yy,phi,[0.0])
            plt.plot(xx, yy, 'k-', lw=0.5, alpha=0.5)
            plt.plot(xx.T, yy.T, 'k-', lw=0.5, alpha=0.5)
            plt.title(pvar)
            plt.pause(.001)



# Curve test.  Write file and print its name at the end
if test:
    v = ss.PyMPI.zbar( ss.variables[pvar].data )
    ny = ss.PyMPI.ny
    v1d =  v[:,int(ny/2)]
    x1d = xx[:,int(ny/2)]
    fname = testName + '.dat'
    numpy.savetxt( fname  , (x1d,v1d) )
    print(fname)
