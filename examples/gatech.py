import numpy 
from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep


## Define a mesh
problem = 'GAtechTest3'             # Problem name: used for data output
nx = 50                             # Points in x
ny = 100                            # Points in y
dim = 2                             # Dimension of the problem
open_angle  = 45.0                  # Total opening angle in degrees
open_angle *= (numpy.pi / 180.0)    # Convert to radians
L_min = 2.62                        # Width of tube at inner radius
h = 40.0                            # Height is y direction
intH = 21.0                         # Height of interface
mwL = 1.0                           # Molecular weight of light
mwH = 3.14                          # Molecular weight of heavy
myGamma = 1.4                       # Gamma of gases (single gamma)

intAmp = .01*h
intFreq = 4.0

# Farfield properties
p0   = 1.0
rho0 = 1.0 #* mwH/mwL

import sys
try:
    res = int( sys.argv[1] )
    nx *= res
    ny *= res
except:
    pass



# Compute max length (x-dir) based on angle/height
L_max = 2.0 * ( L_min/2.0 + h*numpy.tan( open_angle/2.0 )   )

# Define a function of the mesh
def convTube(i,j,k):
    y   = float(j)*h/float(ny-1)
    wgt = float(ny-1-j)/float(ny-1)
    Lh  = L_min*wgt + L_max*(1.0-wgt)
    x   = -Lh/2.0 + Lh*float(i)/(nx-1)
    z = 0
    return x,y,z

dy0 = L_min / float(nx)
jo = (ny-1)*(L_min/L_max)
alpha = h / ( (ny-1 +jo)**2 - jo**2 )
#alpha = h / ( ((ny-1)+jo)**2 - jo**2 )
beta  = - alpha * jo**2
def adapt(i,j,k):
    y = alpha*(j+jo)**2 + beta
    wgt = y/h
    Lh  = L_min*(1.0-wgt) + L_max*wgt
    x   = -Lh/2.0 + Lh*float(i)/(nx-1)
    z = 0
    return x,y,z

# Define a python dictionary for the mesh/problem
mesh_options = {}
mesh_options['coordsys'] = 3          # This assumes curvilinear grids
mesh_options['function'] = adapt      # Function takes (i,j,k) args and returns [x,y,z]
mesh_options['periodicGrid'] = False  # Is the grid itself periodic?  
mesh_options['dim'] = dim             # Dimension of the problem
mesh_options['nn'] = [ nx, ny , 1 ]   # Grid spacing
mesh_options['periodic'] = [False, False, False]

# Initialize a simulation object on a mesh
ss = pyrandaSim(problem,mesh_options)

# Add packages to simulation object
pyBC = pyrandaBC(ss)
ss.addPackage( pyBC )      # BC package allows for "bc.*" functions
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
:Ysgs:      =  1.0e2*(abs(:Yh:) - 1.0 + abs(1.0-:Yh: ) )*gridLen**2
:adiff:     =  gbar( :rho:*numpy.maximum(:Dsgs:,:Ysgs:) / :dt: )
[:Yx:,:Yy:,:Yz:] = grad( :Yh: )
:Jx:        =  :adiff:*:Yx:
:Jy:        =  :adiff:*:Yy:
# Apply boundary conditions
bc.extrap(['rho','p'],['x1','xn','y1'] ,order=1)
bc.extrap(['Yh'],['x1','xn','y1','yn'], order=1 )
bc.slip( [ ['u','v'] ] , ['x1','xn','y1'] )
# Sponge outflow
#bc.const(['p'],['yn'],p0)
#bc.const(['rho'],['yn'],rho0)
#bc.const(['u','v'],['yn'],0.0)
:wgt: = ( 1.0 + tanh( (meshy-myH*(1.0-0.025))/ (.025*myH) ) ) * 0.5
:u: = :u:*(1-:wgt:) + gbary(:u:)*:wgt:
:v: = :v:*(1-:wgt:) + gbary(:v:)*:wgt:
:rho: = :rho:*(1-:wgt:) + gbary(:rho:)*:wgt:
:p:   = :p:  *(1-:wgt:) + gbary(:p:)  *:wgt:
#:w: = :w:*(1-:wgt:) + gbary(:w:)*:wgt:
bc.extrap(['u','v','p','rho'],['yn'])
#bc.extrap( ['u'] , ['yn'], 0.0 )
#bc.exit( ['u'] , ['yn'] )
#bc.exit( ['u'] , ['yn'] , True )
#:vout: = numpy.maximum(0.0,:v:)
#bc.field( ['v'] , ['yn'], :vout: )
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:Yl:  = 1.0 - :Yh:
:rhoYh: = :rho:*:Yh:
:rhoYl: = :rho:*:Yl:
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
# Compute some max time steps
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dtC: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.2* dt.diff(:beta:,:rho:)
:dt:  = numpy.minimum(:dtC:,:dtB:)
:dtM: = 0.2 * dt.diff(:mu:,:rho:)
:dt: = numpy.minimum(:dt:,:dtM:)
:dtY: = 0.2 * dt.diff(:adiff:,:rho:)
:dt: = numpy.minimum(:dt:,:dtY:)
:umag: = sqrt( :u:*:u: + :v:*:v: )
"""

# Add the EOM to the solver
eomDict = {'mwH':mwH,'mwL':mwL,'myH':h}
eomDict['p0'] = p0
eomDict['rho0'] = rho0
ss.EOM(eom, eomDict )


# Initial conditions of the flow
ic = """
# Interaface perturbations and smootings
sinW = height + AMP*sin( (meshx+intWidth/2.0) / intWidth * 2.0 * pi * FREQ )
:Yh: = 0.5 * (1.0 - tanh( (meshy - sinW ) / thick ) )
:Yl: = 1.0 - :Yh:
# Mixed EOS (partial pressures method)
:gamma: = myGamma
:mw: = 1.0 / ( :Yh: / mwH + :Yl: / mwL )
:R: = 1.0 / :mw:
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
# Energy pill for detonation wave
rad = sqrt( (meshx-0.0)**2  +  (meshy-0.0)**2 ) 
:p:  =  1.0 + 200.0 * exp( - rad**2 / 0.5**2 )
:T:  = :p:  #**( (:gamma:-1.0) / :gamma: )
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
intR = intH / h
icDict['height']  = h * intR
icDict['intWidth']= L_min*intR + L_max*(1.0-intR)
icDict['thick'] = h / ny * 2
icDict['AMP'] =  intAmp
icDict['FREQ'] = intFreq

# Set the initial conditions
ss.setIC(ic,icDict)
    
# Write a time loop
time = 0.0
viz = True

# Stopping time of simulation
tt = 200.0

# Variables to write to viz.
wvars = ['Yh','rho','u','v','p','beta','kappa','adiff','mu']

# Start time loop
dump_freq = 3000
CFL = 1.0
dt = ss.variables['dt'].data * CFL*.01

viz_freq = tt / 200.0
viz_dump = viz_freq

ss.plot.figure(1)
ss.plot.clf()
ss.plot.contourf('p',64 )
ss.write( wvars )

while time < tt :
    
    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = min( ss.variables['dt'].data * CFL, dt*1.01)
    dt = min(dt, (tt - time) )

    stab_type = ''
    # Simulation heart-beat
    if dt == ss.variables['dtB'].data:
        stab_type = 'bulk'
    elif dt == ss.variables['dtY'].data:
        stab_type = 'diffusion'
    elif dt == ss.variables['dtM'].data:
        stab_type = 'shear'
    elif dt == ss.variables['dtC'].data:
        stab_type = "CFL"
    else:
        stab_type = "ramp"
        
    ss.iprint("Cycle: %5d --- Time: %10.4e --- deltat: %10.4e ( %s )" % (ss.cycle,time,dt,stab_type)  )
    
    # Simulation heart-beat
    ss.iprint("Cycle: %5d --- Time: %10.4e --- deltat: %10.4e" % (ss.cycle,time,dt)  )


    # Constant time
    if time > viz_dump:
        ss.write( wvars )
        viz_dump += viz_freq

        ss.plot.figure(1)
        ss.plot.clf()
        ss.plot.contourf('p',64 )    

        
    if (ss.cycle%dump_freq == 0) :
        ss.writeRestart()

