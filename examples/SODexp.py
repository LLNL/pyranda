import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC
from pyranda.pyranda import pyrandaRestart


## Define a mesh
L = numpy.pi * 2.0
Npts = 200
Lp = L * (Npts-1.0) / Npts

imesh = """
xdom = (0.0, Lp, Npts)
""".replace('Lp',str(Lp)).replace('Npts',str(Npts))

# Initialize a simulation object on a mesh
ssE = pyrandaSim('sod',imesh)
ssE.addPackage( pyrandaBC(ssE) )

ssC = pyrandaSim('sod',imesh)
ssC.addPackage( pyrandaBC(ssC) )


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -ddx8e(:rho:*:u:)
ddt(:rhou:) =  -ddx8e(:rhou:*:u: + :p: - :tau:)
ddt(:Et:)   =  -ddx8e( (:Et: + :p: - :tau:)*:u: - :qx: )
# Conservative filter of the EoM
#:rho:       =  fbar6e( :rho:  )
#:rhou:      =  fbar6e( :rhou: )
#:Et:        =  fbar6e( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
:T: = :p: / :rho:
# Artificial bulk viscosity (old school way)
#:div:       =  ddx6e(:u:) 
:div:       =  ddx2e(:u:) 
:beta:      =  gbar( abs(dd4x(:div:)*gridLen**2) * :rho:) * CB
:tau:       =  :beta:*:div:
:tx:        =  ddx6e(:T:)
:e:         =  :p: / (:gamma: - 1.0 ) / :rho:
:kappa:     =  gbar( abs(dd4x(:e:)) * :rho: /:T: * :cs: * gridLen )*CK
:qx:        =  :kappa:*:tx:
# Apply constant BCs
bc.extrap(['rho','Et'],['x1'])
bc.const(['u'],['x1','xn'],0.0)
:cs: = sqrt( :p: / :rho: / :gamma: )
""".replace("gridLen",str(Lp/Npts))



#eomEXP = eom.replace('CB',str(10.0)).replace('CK',str(0.5))
#eomCOM = eom.replace('8e','6e').replace('CB',str(1.0)).replace('CK',str(0.05))

eomEXP = eom.replace('8e','6e').replace('CB',str(10.0)).replace('CK',str(0.5))
eomCOM = eom.replace('8e','').replace('CB',str(1.0)).replace('CK',str(0.05))

# Add the EOM to the solver
ssE.EOM(eomEXP)
ssC.EOM(eomCOM)


# Initial conditions SOD shock tube in 1d
ic = """
:gamma: = 1.4
:Et:  = gbar( where( meshx < pi, 1.0/(:gamma:-1.0) , .1 /(:gamma:-1.0) ) )
:rho: = gbar( where( meshx < pi, 1.0    , .125 ) )
"""

# Set the initial conditions
ssE.setIC(ic)
ssC.setIC(ic)
    

# Write a time loop
time = 0.0

# Approx a max dt and stopping time
v = 1.0
dt_max = v / ssE.mesh.nn[0] * 0.75
tt = L/v * .25 #dt_max

# Start time loop
dt = dt_max
cnt = 1
viz_freq = 25
pvar = 'rho'
viz = True

while tt > time:

    # Update the EOM and get next dt
    timeDUM = ssE.rk4(time,dt)
    time    = ssC.rk4(time,dt)


    if cnt%1 == 0:
        ssE.parse(":rho:       =  fbar8e( :rho:  )")
        ssE.parse(":rhou:      =  fbar8e( :rhou: )")
        ssE.parse(":Et:        =  fbar8e( :Et:   )")
        ssE.updateVars()

        ssC.parse(":rho:       =  fbar( :rho:  )")
        ssC.parse(":rhou:      =  fbar( :rhou: )")
        ssC.parse(":Et:        =  fbar( :Et:   )")
        ssC.updateVars()

    
    dt = min(dt_max, (tt - time) )
    
    # Print some output
    ssC.iprint("%s -- %s" % (cnt,time)  )
    cnt += 1
    if viz:

        if (cnt%viz_freq == 0):
            ssE.plot.figure(1)
            plt.clf()
            ssE.plot.plot(pvar,'b.-')
            ssC.plot.plot(pvar,'k.-')

        
#ss.writeGrid()
#ss.write()
