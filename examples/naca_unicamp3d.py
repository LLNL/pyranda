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
    Npts = 100

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

problem = 'naca_test3'

Lp = L * (Npts-1.0) / Npts


# Make and load an NACA mesh
from omesh import naca_omesh
NACA = '2412'
nx = 600
ny = 100
[xnaca,ynaca] = naca_omesh(NACA,nx,ny,
                           PLOT=False,
                           dr0=.0025,fact=1.05,iS=5,
                           te=.08,teSig=12,dratio=10)


dz = .01 

def nacaMesh(i,j,k):
    x = xnaca[i,j]
    y = ynaca[i,j]
    #z = 0.0
    z = dz * float(k)
    return x,y,z

mesh_options = {}
mesh_options['coordsys'] = 3
mesh_options['function'] = nacaMesh
mesh_options['periodic'] = numpy.array([True, False, True])
mesh_options['periodicGrid'] = True  # really should be NOT periodic mesh
mesh_options['dim'] = 3
mesh_options['x1'] = [ -2*Lp , -2*Lp  ,  0.0 ]
mesh_options['xn'] = [ 2*Lp   , 2*Lp    ,  Lp ]
mesh_options['nn'] = [ nx, ny ,  16  ]


# Initialize a simulation object on a mesh
ss = pyrandaSim(problem,mesh_options)
ss.addPackage( pyrandaBC(ss) )
ss.addPackage( pyrandaIBM(ss) )
ss.addPackage( pyrandaTimestep(ss) )


rho0 = 1.0
p0   = 1.0
gamma = 1.4
mach = 0.7
s0 = numpy.sqrt( p0 / rho0 * gamma )
u0 = s0 * mach
e0 = p0/(gamma-1.0) + rho0*.5*u0*u0

Re = 10000
mu0 = u0 * rho0 * 1.0 / Re

# Define the equations of motion
from equation_library import euler_3d
eom = euler_3d
eom +="""
# Apply constant BCs
bc.extrap(['u','v','rho','p'],['yn'])
bc.const(['u'],['yn'],#u0#)
bc.extrap(['rho','p'],['y1'],order=1)
bc.const( ['u','v']  ,['y1'],0.0)
# Sponge outflow
:wgt: = ( 1.0 + tanh( (:rad:- 5.8 )/ (.25) ) ) * 0.5
:u: = :u:*(1-:wgt:) + gbar(:u:)*:wgt:
:v: = :v:*(1-:wgt:) + gbar(:v:)*:wgt:
:w: = :w:*(1-:wgt:) + gbar(:w:)*:wgt:
bc.extrap(['u','v','w'],['xn'])
# Update the conserved quantities
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dtC: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.2* dt.diff(:beta:,:rho:)
:dt: = numpy.minimum(:dtC:,:dtB:)
:dtM: = 0.2* dt.diff(:mu:,:rho:)
:dt: = numpy.minimum(:dt:,:dtM:)
:umag: = sqrt( :u:*:u: + :v:*:v: + :w:*:w:)
"""
eom = eom.replace('#u0#',str(u0)).replace('mu0',str(mu0))


# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = """
:gamma: = 1.4
:R: = 1.0
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
:rho: = 1.0 + 3d()
:p:  =  1.0 + 3d()
:rad: = sqrt( meshx**2 + meshy**2 ) 
:u: = mach * sqrt( :p: / :rho: * :gamma:)
bc.const( ['u'] ,['y1'], 0.0)
:u: = gbar(gbar(gbar(gbar(gbar(gbar( gbar( :u: ) ) )))))
:v: = gbar(gbar(gbar(gbar(gbar(gbar( gbar( :v: ) ) )))))
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""
ic = ic.replace('mach',str(mach))

# Set the initial conditions
ss.setIC(ic)
    

# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
tt = 18.0 #


# Start time loop
cnt = 1
viz_freq = 100
pvar = 'umag'


CFL = 0.5
dt = ss.variables['dt'].data * CFL



v = ss.PyMPI.zbar( ss.variables[pvar].data )
p = ss.PyMPI.zbar( ss.variables['p'].data )

if (not test) :

    ss.plot.figure(2)
    ss.plot.clf()            
    ss.plot.contourf( pvar, 64 , cmap=cm.jet)
    ss.plot.showGrid()    
    

wvars = ['p','rho','u','v','w','T','Et','beta','mu','umag']

ss.write(wvars)

#import pdb
#pdb.set_trace()

ramp = .01
dt *= ramp
while tt > time:
    
    # Update the EOM and get next dt


    
    time = ss.rk4(time,dt)
    dt = ss.variables['dt'].data * CFL
    dt = min(dt, (tt - time) )

    dt *= ramp
    ramp = min( 1.01*ramp, 1.0)
    #ss.iprint(ramp)
    
    

    
    
    # airfoil extrap
    tramp1 = .01
    tramp2 = .05
    if time < tramp1:
        ss.parse(':u: = gbar(:u:)')
        ss.parse(':v: = gbar(:v:)')
        ss.parse(':rho: = gbar(:rho:)')
        ss.parse(':rhou: = :rho:*:u:')
        ss.parse(':rhov: = :rho:*:v:')


    if time > tramp1 and time < tramp2:
        wgt = 1.0 - (time - tramp1) / (tramp2-tramp1)  # Linear ramp
        ss.parse(':u:   = :u:*(1.0-%s)   + %s*gbar(:u:)'   % (wgt,wgt) )
        ss.parse(':v:   = :v:*(1.0-%s)   + %s*gbar(:v:)'   % (wgt,wgt))
        ss.parse(':rho: = :rho:*(1.0-%s) + %s*gbar(:rho:)' % (wgt,wgt))
        ss.parse(':rhou: = :rho:*:u:')
        ss.parse(':rhov: = :rho:*:v:')

    
    
    # Print some output
    stab = 'None'
    if ( ss.variables['dt'].data == ss.variables['dtB'].data):
        stab = "bulk"
    if ( ss.variables['dt'].data == ss.variables['dtM'].data):
        stab = "shear"
    if ( ss.variables['dt'].data == ss.variables['dtC'].data):
        stab = "CFL"
                    
    ss.iprint("%s -- %f --- %f --- %s" % (cnt,time,dt,stab)  )
    cnt += 1
    if viz and (not test):
        v = ss.PyMPI.zbar( ss.variables[pvar].data )
        if (cnt%viz_freq == 1) :#or True:

            #ss.plot.figure(2)
            #ss.plot.clf()            
            #ss.plot.contourf( pvar ,64 , cmap=cm.jet)
            #if ss.PyMPI.master:
            #    plt.xlim([-1,2])
            #    plt.ylim([-2,2])
            #    plt.pause(.01)
                
            #ss.plot.showGrid()  

            ss.write(wvars)



