from __future__ import print_function
import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep
from pyranda.pyrandaMesh import defaultMeshOptions


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
dx2 = L / (2.*Npts - 1.0)
gamma = 1.4
problem = 'sod2dRZ'
mesh_options = defaultMeshOptions()
mesh_options['coordsys'] = 0
mesh_options['periodic'] = numpy.array([False, False, True])
#mesh_options['dim'] = 2
mesh_options['x1'] = [ dx2 , -0.1  ,  0.0 ]
mesh_options['xn'] = [ L   ,  0.1  ,  L ]
mesh_options['nn'] = [ Npts, 1 ,  Npts  ]
mesh_options['symmetric'][0][0] = True
mesh_options['symmetric'][0][1] = True

# Initialize a simulation object on a mesh
ss = pyrandaSim('advection',mesh_options)
ss.addPackage( pyrandaBC(ss) )
ss.addPackage( pyrandaTimestep(ss) )

# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -:mass:
ddt(:rhou:) =  -:xmom:
ddt(:rhow:) =  -:zmom:
ddt(:Et:)   =  -div( (:Et: + :p: - :tau:)*:u: - :tx:*:kappa: ,:v:, (:Et: + :p: - :tau:)*:w: - :tz:*:kappa: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou:, comp=1 )   # Adding in COMP currently breaks this case.
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:w:         =  :rhow: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :w:*:w:) ) * ( :gamma: - 1.0 )
:T:         = :p: / (:rho: * :R: )
:cs:  = sqrt( abs(:p: / :rho: * :gamma: ) )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:,:v:,:w:)
:beta:      =  gbar(abs(ring(:div:))) * :rho: * 7.0e-2
:tau:       =  :beta:*:div:
[:tx:,:ty:,:tz:] = grad(:T:)
:IT:   = 1.0 / :T:
:kappa:     = gbar( :rho: * :cs:**3 * ring(:IT:) ) * 1.0e-3
# Compute divT of momentum
:mass: = div(:rho:*:u:               , :v:, :rho:*:w:)
[:fxx:,:fxy:,:fxz:] = [:rhou:*:u: + :p: - :tau:, :v:, :rhou:*:w:]
[:fyx:,:fyy:,:fyz:] = [      :v:               , :p:-:tau: , :v:       ]
[:fzx:,:fzy:,:fzz:] = [:rhow:*:u:,               :v:, :rhow:*:w: + :p: - :tau:]
[:xmom:,:ymom:,:zmom:] = divT(:fxx:,:fxy:,:fxz:,:fyx:,:fyy:,:fyz:,:fzx:,:fzy:,:fzz:)
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.2* dt.diff(:beta:,:rho:)
:dtT: = 0.2* dt.diff(:kappa:,:rho:*:cv:/:T:)
:dt: = numpy.minimum(:dt:,:dtB:)
:dt: = numpy.minimum(:dt:,:dtT:)
:totMass: = sum( :rho:*cellVol )
:totRmom: = sum( :rho:*:u:*cellVol )
:totZmom: = sum( :rho:*:w:*cellVol )
:totE:    = sum( :Et:*cellVol )
"""

eom += """# Apply constant BCs
#bc.extrap(['rho','Et'],['xn'],order=1)
#bc.const(['u','w'],['xn'],0.0)

:rhou: = :rho:*:u:
:rhow: = :rho:*:w:

"""

print(eom)

# Add the EOM to the solver
ss.EOM(eom)


# Initialize variables
ic = "rad = sqrt( (meshx-pi/2)**2  +  (meshz-pi)**2 ) "


# SOD shock tube in 1d and 2d
ic += """
:gamma: = 1.4
wgt = .5*(1.0-tanh( (rad-pi/4.0)/.1) )   # [0,1]
:Et:  = 1.0/(:gamma:-1.0) * (wgt*.9 + e2)
#:Et:  = gbar( where( rad < pi/2.0, 1.0/(:gamma:-1.0) , .1 /(:gamma:-1.0) ) )
#:rho: = gbar( where( rad < pi/2.0, 1.0    , .125 ) )
:rho: = 1.0*wgt + (1.0-wgt)*rho2

:R: = 1.0
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:

:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :w:*:w:) ) * ( :gamma: - 1.0 )

:cs:  = sqrt( abs(:p: / :rho: * :gamma: ) )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)

"""
rho2 = .125
e2   = .1
parms = {}
parms['rho2'] = rho2
parms['e2'] = e2


# Set the initial conditions
ss.setIC(ic,icDict=parms)
    
# Write a time loop
time = 0.0
viz = True

# Approx a max dt and stopping time
v = 1.0
dt = ss.variables['dt'].data * .1
tt = 2.0

if test:
    tt = 0.5
    
# Start time loop
cnt = 1
viz_freq = 10


cs0 = ss.variables['cs'].data.max()

mass0 = ss.variables['totMass'].data
rmom0 = ss.variables['totRmom'].data + rho2*cs0**2 
zmom0 = ss.variables['totZmom'].data + rho2*cs0**2 
Et0   = ss.variables['totE'].data

pvar = 'T'

while tt > time:

    # Update the EOM and get next dt
    time = ss.rk4(time,dt)
    dt = min(ss.variables['dt'].data, 1.1*dt)
    dt = min(dt, (tt - time) )

    # Print some output
    mass = (mass0 - ss.variables['totMass'].data)/mass0
    rmom = (rmom0 - ss.variables['totRmom'].data)/rmom0
    zmom = (zmom0 - ss.variables['totZmom'].data)/zmom0
    Et   = (Et0   - ss.variables['totE'].data)   /Et0
    
    ss.iprint("%s -- %s   ---  MASS:  %5.4e ---  Rmom:  %5.4e ---  Zmom:  %5.4e ---  Energy:  %5.4e" % (cnt,time,mass,rmom,zmom,Et )  )
    cnt += 1
    if viz and (not test):
        if (cnt%viz_freq == 0) and True:
            ss.plot.figure(1)
            ss.plot.clf()
            ss.plot.plot(pvar,slice2d="j=0;k=50",style='k.-')
            ss.plot.figure(2)
            ss.plot.clf()            
            ss.plot.contourf(pvar,64)



if test:
    half = int(Npts/2.0)
    x = ss.mesh.coords[0][:,0,0]
    p = ss.var('p')[half,0,:]
    fname = testName + '.dat'
    numpy.savetxt( fname  , (x,p) )
    print(fname)
    
