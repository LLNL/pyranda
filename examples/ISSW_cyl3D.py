from __future__ import print_function
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM 
from pyranda.pyranda import pyrandaRestart


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except:
    Npts = 64

try:
    restart = bool(int(sys.argv[2]))
except:
    restart = False


try:
    dim = int(sys.argv[3])
except:
    dim = 3
    

## Define a mesh
L = numpy.pi * 2.0  
gamma = 1.4

problem = 'cyl_ISSW_fine2d'

Lp = L * (Npts-1.0) / Npts
dx_ish = 4.*Lp / float(Npts)
Nz = 72
Lz = Nz * dx_ish

mesh_options = {}
#mesh_options['coordsys'] = 3
#mesh_options['function'] = zoomMesh
mesh_options['periodic'] = numpy.array([False, False, True])
mesh_options['dim'] = 3
mesh_options['x1'] = [ -2*Lp + 6., -2*Lp  ,  0.0 ]
mesh_options['xn'] = [ 2*Lp  + 6. , 2*Lp    ,  Lz ]
mesh_options['nn'] = [ Npts, Npts ,  Nz  ]
if dim == 2:
    mesh_options['nn'] = [ Npts, Npts ,  1  ]

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
ddt(:rho:)  =  -div(:rho:*:u:,  :rho:*:v:, :rho:*:w:)
ddt(:rhou:) =  -div(:rhou:*:u: + :p: - :tau:, :rhou:*:v:, :rhou:*:w:)
ddt(:rhov:) =  -div(:rhov:*:u:, :rhov:*:v: + :p: - :tau:, :rhov:*:w:)
ddt(:rhow:) =  -div(:rhow:*:u:, :rhow:*:v:, :rhow:*:w: + :p: - :tau:)
ddt(:Et:)   =  -div( (:Et: + :p: - :tau:)*:u: - :tx:*:kappa: ,(:Et: + :p: - :tau:)*:v: - :ty:*:kappa: , (:Et: + :p: - :tau:)*:w: - :tz:*:kappa: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:w:         =  :rhow: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:) ) * ( :gamma: - 1.0 )
:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:,:v:,:w:)
:beta:      =  gbar( ring(:div:) * :rho: ) * 7.0e-3
:tau:       =  :beta: * :div: 
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
# Apply constant BCs
[:u:,:v:,:w:] = ibmV( [:u:,:v:,:w:], :phi:, [:gx:,:gy:,:gz:] )
:rho: = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
:p:   = ibmS( :p:   , :phi:, [:gx:,:gy:,:gz:] )
bc.extrap(['rho','p','u'],['xn'])
bc.const(['u'],['x1','y1','yn'],u0)
bc.const(['v'],['x1','xn','y1','yn'],0.0)
bc.const(['w'],['x1','xn','y1','yn'],0.0)
bc.const(['rho'],['x1','y1','yn'],rho0)
bc.const(['p'],['x1','y1','yn'],p0)
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*.6
:dtB: = 0.05* dt.diff(:beta:,:rho:)
:dt: = numpy.minimum(:dt:,:dtB:)
:umag: = sqrt( :u:*:u: + :v:*:v: + :w:*:w: )
"""
eom = eom.replace('u0',str(u0)).replace('p0',str(p0)).replace('rho0',str(rho0))


# Initialize variables
ic = """
:gamma: = 1.4
:R: = 1.0
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
#rad = sqrt( (meshx-pi)**2  +  (meshy-pi)**2 ) 
rad = sqrt( meshx**2  +  meshy**2 ) 
:phi: = rad - pi/4.0
rad = sqrt( (meshx+3)**2  +  (meshy-5.0)**2 ) 
:phi: = numpy.minimum(rad - pi/4.0,:phi:)
rad = sqrt( (meshx-3)**2  +  (meshy+5.0)**2 ) 
:phi: = numpy.minimum(rad - pi/4.0,:phi:)
:rho: = 1.0 + 3d()
:p:  =  1.0 + 3d() #exp( -(meshx-1.5)**2/.25**2)*.1
:u: = where( :phi:>0.5, mach * sqrt( :p: / :rho: * :gamma:) , 0.0 )
#:u: = mach * sqrt( :p: / :rho: * :gamma:)
:u: = gbar( gbar( :u: ) )
:v: = 0.0 + 3d()
:w: = 0.0 + 3d()
:w: = (numpy.random.random( :w:.shape ) - 0.5)*.25
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
[:gx:,:gy:,:gz:] = grad( :phi: )
:gx: = gbar( :gx: )
:gy: = gbar( :gy: )
"""
ic = ic.replace('mach',str(mach))


# Option to pick up a restart
CFL = 1.0
outVars = ['p','u','v','w','phi','rho','T']

# Approx a max dt and stopping time
tt = 25.0 #

# Start time loop
cnt = 1

# Viz loop
viz_N = 500
viz_times = numpy.linspace(0,tt,viz_N)
viz = False  # Interactive plots?

# Restart
restart_freq = 10


# Pyranda Simulation
if restart:
    ss = pyrandaRestart(problem)
    dt = ss.deltat

else:
    # Initialize a simulation object on a mesh
    ss = pyrandaSim(problem,mesh_options)
    ss.addPackage( pyrandaBC(ss) )
    ss.addPackage( pyrandaIBM(ss) )
    ss.addPackage( pyrandaTimestep(ss) )

    # Add the EOM/IC to the solver
    ss.EOM(eom)
    ss.setIC(ic)
    dt = ss.variables['dt'].data * CFL*.01

    ss.write(outVars)
    

time = ss.time
while tt > time:
    
    # Update the EOM and get next dt
    t_old = time*1.0
    time = ss.rk4(time,dt)
    dt = min( ss.variables['dt'].data * CFL, dt*1.1)
    dt = min(dt, (tt - time) )

    
    # Print some output
    cnt = ss.cycle
    ss.iprint("%s -- %s --- %f" % (cnt,time,dt)  )

    # Write dumps
    if (cnt%restart_freq == 0):
        ss.iprint("Restart file written at cycle %s" % cnt)
        ss.writeRestart()

    # Write/plot viz
    for vtime in viz_times:
        if (t_old < vtime) and (time >= vtime):
            ss.iprint("Viz file written at t=%s" % time)
            ss.write(outVars)
            if viz :
                ss.plot.clf()            
                ss.plot.contourf('p',64,cmap=cm.jet)



            

