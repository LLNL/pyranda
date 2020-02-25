import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaIBM

import matplotlib


# Issues... Skeleton point at center of level set. Forcing function cant get to unit grad vector


def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)
    #plt.show()


def plotFix( plots ):
    px_max = 1300
    py_max = 1000
    """
    Arrange the plots in reasonable array
    """
    xw = 640
    yw = 550
    x1 = y1 = 0
    for p in plots:
        #print(x1)
        #print(y1)
        move_figure(p,x1,y1)
        x1 += xw
        if x1 > px_max:
            x1 = 0
            y1 += yw
            if y1 > py_max:
                y1 = 0
    
## Define a mesh
L = numpy.pi * 2.0
Npts = 100
Lp = L * (Npts-1.0) / Npts

imesh = """
xdom = (0.0, Lp, Npts)
ydom = (0.0, Lp, Npts)
""".replace('Lp',str(Lp)).replace('Npts',str(Npts))

# Initialize a simulation object on a mesh
ss = pyrandaSim('bubble',imesh)
ss.addPackage( pyrandaBC(ss) )
ss.addPackage( pyrandaIBM(ss) )

deltaPhi = L / float(Npts)
# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rhoA:)  =  -ddx(:rhoA:*:uA:) - ddy(:rhoA:*:vA:)
ddt(:rhoB:)  =  -ddx(:rhoB:*:uB:) - ddy(:rhoB:*:vB:)
ddt(:rhouA:) =  -ddx(:rhouA:*:uA: + :pA: - :tauA:) -ddy(:rhouA:*:vA:) - :FAx:
ddt(:rhovA:) =  -ddx(:rhovA:*:uA:) - ddy(:rhovA:*:vA: + :pA: - :tauA:) - :FAy:
ddt(:rhouB:) =  -ddx(:rhouB:*:uB: + :pB: - :tauB:)-ddy(:rhouB:*:vB:) - :FBx: 
ddt(:rhovB:) =  -ddx(:rhovB:*:uB:) - ddy(:rhovB:*:vB: + :pB: - :tauB:) - :FBy:
ddt(:EtA:)   =  -ddx( (:EtA: + :pA: - :tauA:)*:uA: ) - ddy( (:EtA: + :pA: - :tauA:)*:vA: ) - :FAx:*:uA: - :FAy:*:vA:
ddt(:EtB:)   =  -ddx( (:EtB: + :pB: - :tauB:)*:uB: ) - ddy( (:EtB: + :pB: - :tauB:)*:vB: ) - :FBx:*:uB: - :FBy:*:vB:
ddt(:phi:)   = -:gx:*:uphi: -:gy:*:vphi: + sign(:phi:) * where(abs(:phi:)>4.*deltaPhi,1.0,0.0)*( .01*(1.0-:mgp:) * deltaPhi / (deltat+1.0e-10))     #+ 1.0e-4 * :artPhi: * :curv: * deltaPhi**2 / (deltat) ) 
#ddt(:uphi:)  =  -ddx( :uphi:*:uphi: - :tauphi:) - :Fphi:/(:rhoA: + :rhoB:)
ddt(:uphi:)  =   - :Fphix:/:rhophi:
ddt(:vphi:)  =   - :Fphiy:/:rhophi:
# Conservative filter of the EoM
:rhoA:       =  fbar( :rhoA:  )
:rhoB:       =  fbar( :rhoB:  )
:rhouA:      =  fbar( :rhouA: )
:rhouB:      =  fbar( :rhouB: )
:rhovA:      =  fbar( :rhovA: )
:rhovB:      =  fbar( :rhovB: )
:EtA:        =  fbar( :EtA:   )
:EtB:        =  fbar( :EtB:   )
# Phi gradient
#:phi: = ibmRD( :phi: )
:phi: = where( :phi: < -2*deltaPhi, gbar(:phi:),:phi:)
[:gx:,:gy:,:gz:] = grad( :phi: )
#:gx: = gbar(:gx:)
#:gy: = gbar(:gy:)
:curv: = abs(div(:gx:,:gy:))  
:mgp:        =   sqrt( :gx:**2 + :gy:**2 )
:uphimag: = sqrt( :uphi:*:uphi: + :vphi:*:vphi: )
#:artPhi: = gbar( abs( ring(:mgp:) ) ) / (deltat+1.0e-10)  # need 1/t left over
:artPhi: = gbar( where( :curv: > .2/deltaPhi , 1.0, 0.0 ) )
:wgt: = .5*(1.0 - tanh( (:phi: - 2.0*deltaPhi)/(deltaPhi) ))
:wgt2: = .5*(1.0 - tanh( (-:phi: - 2.0*deltaPhi)/(deltaPhi) ))
:wgt3: = 1.0 - numpy.minimum( :wgt:, :wgt2: )
#:artPhi: = :artPhi: #*:wgt3:
#:tphi: = .01*gbar(:phi:)   +    :phi:*.99
#:phi: = :tphi:*:wgt3:   +    :phi:*(1.0 - :wgt3:)
# Update the primatives and enforce the EOS
:uA:         =  :rhouA: / :rhoA:
:uB:         =  :rhouB: / :rhoB:
:vA:         =  :rhovA: / :rhoA:
:vB:         =  :rhovB: / :rhoB:
:pA:         =  ( :EtA: - .5*:rhoA:*(:uA:*:uA: + :vA:*:vA:) ) * ( :gamma: - 1.0 )
:pB:         =  ( :EtB: - .5*:rhoB:*(:uB:*:uB: + :vB:*:vB:) ) * ( :gamma: - 1.0 )
:wgt:        = (1.0-tanh(:phi:/deltaPhi)**2)
:dphi:       = :wgt:*.5/ deltaPhi
:Fphix:      = 1.0*:dphi:*:gx:*(:pA:-:pB:)
:Fphiy:      = 1.0*:dphi:*:gy:*(:pA:-:pB:)
:wgt2:       = .5*(1.0 - tanh( :phi:/deltaPhi ))
:rhophi:     = :rhoA: + :wgt2:*(:rhoB:-:rhoA:)
:rhophi:         = ibmS( :rhophi: , :phi:, [:gx:,:gy:,:gz:] )
:rhophi:         = ibmS( :rhophi: , -:phi:, [-:gx:,-:gy:,-:gz:] )
# Close level set velocity
:uphi:         = ibmS( :uphi: , :phi:, [:gx:,:gy:,:gz:] )
:uphi:         = ibmS( :uphi: , -:phi:, [-:gx:,-:gy:,-:gz:] )
# Apply constant BCs
bc.extrap(['rhoA','EtA'],['x1'])
bc.extrap(['rhoB','EtB'],['xn'])
bc.const(['uA'],['x1'],0.0)
bc.const(['uB'],['xn'],0.0)
# IBM
:rhoA:         = ibmS( :rhoA: , :phi:, [:gx:,:gy:,:gz:] )
[:uA:,:vA:,:w:] = ibmV( [:uA:,:vA:,0.0], :phi:, [:gx:,:gy:,:gz:], [:uphi:,:vphi:,0.0] )
:pA:          = ibmS( :pA: , :phi:, [:gx:,:gy:,:gz:] )
:EtA:         = :pA: / ( :gamma: - 1.0 ) + .5*:rhoA:*(:uA:*:uA:+ :vA:*:vA:) 
:rhouA:        = :rhoA:*:uA:
:rhovA:        = :rhoA:*:vA:
:rhoB:         = ibmS( :rhoB: , -:phi:, [-:gx:,-:gy:,-:gz:] )
[:uB:,:vB:,:w:] = ibmV( [:uB:,:vB:,0.0], -:phi:, [-:gx:,-:gy:,-:gz:], [:uphi:,:vphi:,0.0] )
:pB:          = ibmS( :pB: , -:phi:, [-:gx:,-:gy:,-:gz:] )
:EtB:         = :pB: / ( :gamma: - 1.0 ) + .5*:rhoB:*(:uB:*:uB:+ :vB:*:vB:) 
:rhouB:        = :rhoB:*:uB:
:rhovB:        = :rhoB:*:vB:
# Artificial bulk viscosity (old school way)
:divA:       =  ddx(:uA:) + ddy(:vA:)
:divB:       =  ddx(:uB:) + ddy(:vB:)
:betaA:      =  gbar( ring(:divA:) * :rhoA:) * 7.0e-2
:betaB:      =  gbar( ring(:divB:) * :rhoB:) * 7.0e-2
:tauA:       =  :betaA:*:divA:
:tauB:       =  :betaB:*:divB:
:u: = where(:phi: < 0.0 , :uB:, :uA: )
:v: = where(:phi: < 0.0 , :vB:, :vA: )
:rho: = where(:phi: < 0.0 , :rhoB:, :rhoA: )
:p: = where(:phi: < 0.0 , :pB:, :pA: )
""".replace('deltaPhi',str(deltaPhi))

# Add the EOM to the solver
ss.EOM(eom)

# Initial conditions SOD shock tube in 1d
ic = """
:gamma: = 1.4 + 3d()
pert    = .0*exp( -abs(meshx-1.0)**2/(.2**2))
:pA:    =  gbar(where(meshx < pi/2., 1.0, 0.1)) * ( 1.0 + pert)
#:rhoA:  = (1.0 + 3d()) * ( 1.0 + pert*(:gamma:-1) )
:rhoA:    =  gbar(where(meshx < pi/2., 1.0, 0.125))
#:pA: += 1.0
#:rhoA: += 1.0
:uA:    = pert
:rhouA: = :rhoA: * :uA:
:EtA:   = :pA: / ( :gamma: - 1.0 ) + .5*:rhoA:*(:uA:*:uA:) 
#:pB: = gbar(where(meshx < pi/2.0, 0.1, 0.1)) * ( 1.0 + pert)
:pB: += 0.1
:EtB:   =  :pB:/(:gamma:-1.0)
##:rhoB:  = 0.125 + 3d() 
:rhoB:  = 0.05 + 3d() 
rad = sqrt( (meshx-pi)**2 + (meshy-pi)**2 )
:phi:   = rad - pi/4.0
"""

# Set the initial conditions
ss.setIC(ic)
    

# Write a time loop
time = 0.0

# Approx a max dt and stopping time
v = 1.0
dt_max = v / ss.mesh.nn[0] * 0.125
tt = L/v * 4.025 #dt_max

# Start time loop
dt = dt_max
cnt = 1
viz_freq = 30
pvar = 'rhoA'
viz = True

import pdb
#pdb.set_trace()

wvars = ['rho','rhoA','rhoB','uA','vA','uB','vB','phi','u','v','p','artPhi','mgp','curv']

ss.write(wvars)

while tt > time:

    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = min(dt_max, (tt - time) )

    
    # Print some output
    ss.iprint("%s -- %s" % (cnt,time)  )
    cnt += 1

    #if (cnt%viz_freq == 0):
    #    poo = raw_input('Paused: any key to continue')



    if viz:
        if (cnt%viz_freq == 0):
            #updatePlots()


            ss.write(wvars)
                
            ss.plot.figure(1)
            ss.plot.clf()
            ss.plot.contourf('rhoA',32, cmap=cm.jet)
            ss.plot.contour('phi',[0.0])

            ss.plot.figure(2)
            ss.plot.clf()
            ss.plot.contourf('rhoB',32, cmap=cm.jet)
            ss.plot.contour('phi',[0.0])

            ss.plot.figure(3)
            ss.plot.clf()
            ss.plot.contourf('u',32, cmap=cm.jet)
            ss.plot.contour('phi',[0.0])



#ss.writeGrid()
#ss.write(wvars)
