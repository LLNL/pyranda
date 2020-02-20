import sys
import time
import numpy 
import matplotlib.pyplot as plt
from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep


Npts = 100

## Define a mesh
L    = 1.0



mesh_options = {}
mesh_options['coordsys'] = 0
mesh_options['periodic'] = numpy.array([False, False, True])
mesh_options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
mesh_options['xn'] = [ L   , L    ,  L ]
mesh_options['nn'] = [ Npts, Npts ,  1  ]



# Initialize a simulation object on a mesh
pysim = pyrandaSim('quantumHydro',mesh_options)
pysim.addPackage( pyrandaBC(pysim) )
pysim.addPackage( pyrandaTimestep(pysim) )


eom ="""
# Primary Equations of motion here
ddt(:n:)  = - div( :nu: , :nv: )
ddt(:nu:) = - div( :u:*:nu: , :v:*:nu: ) - :n:/:m:*ddx(:phi:+:Q:) - 1.0/:m:*ddx(:Pe:) 
ddt(:nv:) = - div( :u:*:nv: , :v:*:nv: ) - :n:/:m:*ddy(:phi:+:Q:) - 1.0/:m:*ddy(:Pe:)
:Pe: = :hbar:**2 * (3*pi**2)**(2./3)*:n:**(5./3) / (5*:m:)
:Q: = - :hbar:**2 / (2*:m:) * lap( sqrt(:n:) ) / sqrt(:n:)  
Delta(:phi:) = 4*pi*(:n:-:nion:)
#:phi: *= 0.0
:n: = fbar(:n:)
:nu: = fbar(:nu:)
:nv: = fbar(:nv:)
:u: = :nu: / :n:
:v: = :nv: / :n:
# Boundary conditions
bc.const(['u'],['x1','xn'],0.0)
bc.const(['v'],['y1','yn'],0.0)
bc.extrap(['n'],['x1','xn','y1','yn'])
#:dt: = dt.courant(:u:,:v:,0.0,0.0)*.1
"""

pysim.EOM(eom)


# Define the initial conditions here
ic = """
:m: = 1.0
:nion: = 1.0
:hbar: = 1.0
rad = sqrt( (meshx-.5)**2 + (meshy-.5)**2 )
thick = .05
:n: += 1.0 + 1.0e-2 * exp( -rad/thick )
:u: += 0.0
:v: += 0.0
:nu: = :n:*:u: 
:nv: = :n:*:v: 
"""
pysim.setIC(ic)



# Variables to write to viz.
wvars = ['n','u','v','Pe','Q','phi']

tFinal = 1.0
time = 0.0
dt = 1.0e-5
cyc = 0
while tFinal > time:
    cyc += 1
    time = pysim.rk4(time,dt)
    
    print("Cycle: %s" % cyc)
    
    if (cyc%50 == 0):
        pysim.plot.figure(1)
        pysim.plot.clf()
        pysim.plot.contourf('n',32)
    
