import numpy 
from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep


## Define a mesh
problem = 'GAtechTest2'              # Problem name: used for data output
nx = 100                            # Points in x
ny = 100                            # Points in y
dim = 2                             # Dimension of the problem
open_angle  = 45.0                  # Total opening angle in degrees
open_angle *= (numpy.pi / 180.0)    # Convert to radians
L_min = 2.0                         # Width of tube at inner radius
h = 10.0                            # Height is y direction
mwL = 1.0                           # Molecular weight of light
mwH = 3.0                           # Molecular weight of heavy
myGamma = 1.4                       # Gamma of gases (single gamma)

# Compute max length (x-dir) based on angle/height
L_max = 2.0 * ( L_min/2.0 + h*numpy.sin( open_angle ) )

# Define a function of the mesh
def convTube(i,j,k):
    y   = float(j)*h/float(ny-1)
    wgt = float(ny-1-j)/float(ny-1)
    Lh  = L_min*wgt + L_max*(1.0-wgt)
    x   = -Lh/2.0 + Lh*float(i)/(nx-1)
    z = 0
    return x,y,z

# Define a python dictionary for the mesh/problem
mesh_options = {}
mesh_options['coordsys'] = 3          # This assumes curvilinear grids
mesh_options['function'] = convTube   # Function takes (i,j,k) args and returns [x,y,z]
mesh_options['periodicGrid'] = False  # Is the grid itself periodic?  
mesh_options['dim'] = dim             # Dimension of the problem
mesh_options['nn'] = [ nx, ny , 1 ]   # Grid spacing
mesh_options['periodic'] = [False, False, False]

# Initialize a simulation object on a mesh
ss = pyrandaSim(problem,mesh_options)

# Add packages to simulation object
ss.addPackage( pyrandaBC(ss) )          # BC package allows for "bc.*" functions
ss.addPackage( pyrandaTimestep(ss) )    # Timestep package allows for "dt.*" functions

# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rhoYh:) =  -div(:rhoYh:*:u: - :Jx: , :rhoYh:*:v: - :Jy:  )
ddt(:rhoYl:) =  -div(:rhoYl:*:u: + :Jx: , :rhoYl:*:v: + :Jy:  )
ddt(:rhou:)  =  -div(:rhou:*:u: - :tauxx:, :rhou:*:v: - :tauxy:)
ddt(:rhov:)  =  -div(:rhov:*:u: - :tauxy:, :rhov:*:v: - :tauyy:)
ddt(:Et:)    =  -div( (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tx:*:kappa:, (:Et: - :tauyy:)*:v: - :tauxy:*:u: - :ty:*:kappa: )
# Conservative filter of the EoM
:rhoYh:     =  fbar( :rhoYh:  )
:rhoYl:     =  fbar( :rhoYl:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:rho:       = :rhoYh: + :rhoYl:
:Yh:        =  :rhoYh: / :rho:
:Yl:        =  :rhoYl: / :rho:
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v:) ) * ( :gamma: - 1.0 )
:mw: = 1.0 / ( :Yh: / mwH + :Yl: / mwL )
:R: = 1.0 / :mw:
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:,:v:)
[:ux:,:uy:,:tz:] = grad( :u: )
[:vx:,:vy:,:tz:] = grad( :v: )
:S:         = sqrt( :ux:*:ux: + :vy:*:vy: + .5*((:uy:+:vx:)**2))
:beta:      =  gbar( ring(:div:) * :rho: ) * 7.0e-2
:mu:        =  gbar( abs(ring(:S:  )) ) * :rho: * 1.0e-3
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauxy:     =  :mu:*(:uy:+:vx:)
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
# Artificial species diffusivities
:Dsgs:      =  ring(:Yh:) * 1.0e-4
:Ysgs:      =  1.0e2*(abs(:Yh:) - 1.0 + abs(1.0-:Yh: ) )*gridLen
:adiff:     =  gbar( :rho:*numpy.maximum(:Dsgs:,:Ysgs:) / :dt: )
[:Yx:,:Yy:,:Yz:] = grad( :Yh: )
:Jx:        =  :adiff:*:Yx:
:Jy:        =  :adiff:*:Yy:
# Apply boundary conditions
bc.extrap(['rho','p'],['x1','xn','y1','yn'] ,order=1)
bc.extrap(['Yh'],['x1','xn','y1','yn'], order=1 )
bc.const( ['u','v'] , ['x1','xn','y1','yn'], 0.0 )
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:Yl:  = 1.0 - :Yh:
:rhoYh: = :rho:*:Yh:
:rhoYl: = :rho:*:Yl:
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
# Compute some max time steps
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.2* dt.diff(:beta:,:rho:)
:dt:  = numpy.minimum(:dt:,:dtB:)
:dtM: = 0.2 * dt.diff(:mu:,:rho:)
:dt: = numpy.minimum(:dt:,:dtM:)
:umag: = sqrt( :u:*:u: + :v:*:v: )
"""

# Add the EOM to the solver
ss.EOM(eom,{'mwH':mwH,'mwL':mwL})


# Initial conditions of the flow
ic = """
# Interaface perturbations and smootings
sinW = 5.0 + .1*sin( meshx * 4.0 )
:Yh: = 0.5 * (1.0 - tanh( (meshy - sinW ) / .3 ) )
:Yl: = 1.0 - :Yh:
# Mixed EOS (partial pressures method)
:gamma: = myGamma
:mw: = 1.0 / ( :Yh: / mwH + :Yl: / mwL )
:R: = 1.0 / :mw:
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
# Energy pill for detonation wave
rad = sqrt( (meshx-0.0)**2  +  (meshy-1.0)**2 ) 
:p:  =  1.0 + 10.0 * exp( - rad**2 / 0.5**2 )
:T: = :p:
:rho: = :p: / ( :R: * :T: ) 
:u: = 3d()
:v: = 3d()
# Form conservatives from primatives
:rhoYh: = :rho:*:Yh:
:rhoYl: = :rho:*:Yl:
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
# Sound speed and initial time step size
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""
icDict = {}
icDict['myGamma'] = myGamma
icDict['mwH']     = mwH
icDict['mwL']     = mwL

# Set the initial conditions
ss.setIC(ic,icDict)
    
# Write a time loop
time = 0.0
viz = True

# Stopping time of simulation
tt = 30.0

# Variables to write to viz.
wvars = ['Yh','rho','u','v','p','beta','kappa','adiff','mu']

# Start time loop
dump_freq = 3000
CFL = 1.0
dt = ss.variables['dt'].data * CFL*.01

viz_freq = tt / 200.0
viz_dump = viz_freq

while time < tt :
    
    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = min( ss.variables['dt'].data * CFL, dt*1.1)
    dt = min(dt, (tt - time) )
    
    # Simulation heart-beat
    ss.iprint("Cycle: %5d --- Time: %10.4e --- deltat: %10.4e" % (ss.cycle,time,dt)  )


    # Constant time
    if time > viz_dump:
        ss.write( wvars )
        viz_dump += viz_freq
    
    if (ss.cycle%dump_freq == 0) :
        ss.writeRestart()

        #ss.plot.figure(1)
        #ss.plot.clf()
        #ss.plot.contourf('p',64 )    
