import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except:
    Npts = 32

try:
    test = bool(int(sys.argv[2]))
except:
    test = False

problem = 'TGvortex'

## Define a mesh
Npts = 32
imesh = """
xdom = (0.0, 2*pi*FF,  Npts, periodic=True)
ydom = (0.0, 2*pi*FF,  Npts, periodic=True)
zdom = (0.0, 2*pi*FF,  Npts, periodic=True)
""".replace('Npts',str(Npts)).replace('pi',str(numpy.pi)).replace('FF',str( float(Npts-1)/Npts ) )

    
# Initialize a simulation object on a mesh
ss = pyrandaSim(problem,imesh)
ss.addPackage( pyrandaTimestep(ss) )


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)                  - ddy(:rho:*:v:)                  - ddz(:rho:*:w:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tauxx:) - ddy(:rhou:*:v: - :tauxy:)       - ddz(:rhou:*:w: - :tauxz:)
ddt(:rhov:) =  -ddx(:rhov:*:u: - :tauxy:)       - ddy(:rhov:*:v: + :p: - :tauyy:) - ddz(:rhov:*:w: - :tauyz:)
ddt(:rhow:) =  -ddx(:rhow:*:u: - :tauxz:)       - ddy(:rhow:*:v: - :tauyz:)       - ddz(:rhow:*:w: + :p: - :tauzz:)
ddt(:Et:)   =  -ddx( (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tauxz:*:w: ) - ddy( (:Et: - :tauyy:)*:v: -:tauxy:*:u: - :tauyz:*:w:) - ddz( (:Et: - :tauzz:)*:w: - :tauxz:*:u: - :tauyz:*:v: )
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
# Artificial bulk viscosity 
:ux:        =  ddx(:u:)
:vy:        =  ddy(:v:)
:wz:        =  ddz(:w:)
:div:       =  :ux: + :vy: + :wz:
# Remaining cross derivatives
:uy:        =  ddy(:u:)
:uz:        =  ddz(:u:)
:vx:        =  ddx(:v:)
:vz:        =  ddz(:v:)
:wy:        =  ddy(:w:)
:wx:        =  ddx(:w:)
:enst:      = sqrt( (:uy:-:vx:)**2 + (:uz: - :wx:)**2 + (:vz:-:wy:)**2 )
:tke:       = :rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:S:         = sqrt( :ux:*:ux: + :vy:*:vy: + :wz:*:wz: + .5*((:uy:+:vx:)**2 + (:uz: + :wx:)**2 + (:vz:+:wy:)**2) )
:mu:        =  gbar( abs(ring(:S:  )) ) * :rho: * 1.0e-4
:beta:      =  gbar( abs(ring(:div:)) * :rho: )  * 7.0e-3
:taudia:    =  (:beta:-2./3.*:mu:) *:div:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia: - :p:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia: - :p:
:tauzz:     =  2.0*:mu:*:wz:   + :taudia: - :p:
:tauxy:     = :mu:*(:uy:+:vx:) + :taudia:
:tauxz:     = :mu:*(:uz:+:wx:) + :taudia:
:tauyz:     = :mu:*(:vz:+:wz:) + :taudia:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*1.0
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:beta:,:rho:))
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:mu:,:rho:))
"""

# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = """
:gamma: = 1.4
u0 = 1.0
p0 = 100.0
rho0 = 1.0
L = 1.0
:u: =  u0*sin(meshx/L)*cos(meshy/L)*cos(meshz/L)
:v: = -u0*cos(meshx/L)*sin(meshy/L)*cos(meshz/L)
:w: = 0.0*:u:
:p:  = p0 + rho0/16.0*( ( cos(2.*meshx/L) + cos(2.*meshy/L) ) * ( cos(2.*meshz/L) + 2.0 ) - 2.0 )
:rho: = rho0 + 0.0*:u:
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
:rhow: = :rho:*:w:
:Et:  = :p: / (:gamma:-1.0) + 0.5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:cs:  = sqrt( :p: / :rho: * :gamma: )
:tke: = :rho:*(:u:*:u: + :v:*:v: + :w:*:w:)
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
"""

# Set the initial conditions
ss.setIC(ic)
    
# Length scale for art. viscosity
# Initialize variables
x = ss.mesh.coords[0]
y = ss.mesh.coords[1]
z = ss.mesh.coords[2]

# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time


# Mesh for viz on master
xx   =  ss.PyMPI.zbar( x ) / Npts
yy   =  ss.PyMPI.zbar( y ) / Npts
#xx   =   x[:,:,16] / Npts
#yy   =   y[:,:,16] / Npts
ny = ss.PyMPI.ny

# Start time loop
CFL = 0.5
dt = ss.variables['dt'].data * CFL

# Viz
cnt = 1
viz_freq = 200
pvar = 'u'

tke0 = ss.var('tke').sum()
enst0 = ss.var('enst').sum()
TKE = []
ENST = []
TIME = []

tstop = 1.0
if test:
    tstop = .1
    
while time < tstop:

    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = ss.variables['dt'].data * CFL

    # Print some output
    tke = ss.var('tke').sum()/tke0
    enst = ss.var('enst').sum()/enst0

    TIME.append(time)
    ENST.append(enst)
    TKE.append(tke)
    
    
    ss.iprint("%s -- %s --- TKE: %s " % (cnt,time,tke)  )
    cnt += 1
    if viz:

        dd = numpy.where( z == 0.0, 1.0, 0.0 )
        v = ss.PyMPI.zbar( dd*ss.variables[pvar].data )
        
        if (ss.PyMPI.master and (cnt%viz_freq == 0)) and True:
            plt.figure(2)
            plt.clf()            
            plt.contourf( xx,yy,v ,64 , cmap=cm.jet)
            plt.title(pvar)
            plt.pause(.001)

if test:
    print enst
            
#if (ss.PyMPI.master):
#    plt.figure()
#    plt.plot(TIME,TKE,'k--')
#    plt.plot(TIME,ENST,'b--')
#    plt.show()            

