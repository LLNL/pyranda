import sys
import time
import numpy 
import matplotlib.pyplot as plt

from pyranda import pyrandaSim

# Try to get args
try:
    Npts = int(sys.argv[1])
except:
    Npts = 100

try:
    test = bool(int(sys.argv[2]))
except:
    test = False

## Define a mesh
L = numpy.pi * 2.0  
Lp = L * (Npts-1.0) / Npts


# Define the domain/mesh
imesh = """
Lp = %s
Npts = %d
xdom = (0.0, Lp,  Npts, periodic=False)
""" % ( Lp, Npts)

# Initialize a simulation object on a mesh
ss = pyrandaSim('advection',imesh)

# Define the equations of motion
ss.EOM(" ddt(:phi:)  =  -:c: * ddx(:phi:) ")


# Initialize variables
ic = """
r   = sqrt( (meshx-pi)**2  )
:phi: = 1.0 + 0.1 * exp( -(r/(pi/4.0))**2 )
:phi2: = 1.0*:phi:
:c:   = 1.0
"""
ss.setIC(ic)

#ss.variables["u"].data += 1.0

x  = ss.mesh.coords[0].data
xx =  ss.PyMPI.zbar( x )

# Time step size
v = 1.0
dt_max = v / ss.mesh.nn[0] * L * .90
tt = L/v * 1.0 

# Test constant ddx
ss.parse(':phi: = meshx')
ss.plot.plot('phi')
ss.parse(':phi2: = ddx(:phi:)' )
ss.plot.plot('phi2')


ss.plot.figure(2)
ss.parse(':phi: = sin( 4.0*meshx )')
ss.plot.plot('phi')
ss.parse(':phi2: = dd4x(:phi:)' )
ss.plot.plot('phi2')


ss.plot.figure(3)
ss.parse(':phi: = sin( 4.0*meshx )')
ss.plot.plot('phi')
ss.parse(':phi2: = fbar(:phi:)' )
ss.plot.plot('phi2','bo')



ss.plot.figure(4)
ss.plot.plot('phi')
ss.parse(':phi2: = gbar(:phi:)' )
ss.plot.plot('phi2','bo')
