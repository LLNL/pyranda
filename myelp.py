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
ydom = (0.0, Lp,  Npts, periodic=False)
""" % ( Lp, Npts)

# Initialize a simulation object on a mesh
ss = pyrandaSim('Poisson',imesh)

# Define the equations of motion
eom = """
:wx: = meshx * ( 2.0 * pi - meshx )
:wy: = meshy * ( 2.0 * pi - meshy )
:phi: =  cos( meshx ) * :wy: + cos( meshy ) * :wx:
Delta(:ff:) = cos( meshx ) * :wy: + cos( meshy ) * :wx:
"""
ss.EOM(eom)


# Initialize variables
ic = """
#r   = sqrt( (meshx-pi)**2 ) #  + (meshy-pi)**2 )
#:phi: = 0.1 * exp( -(r/(pi/4.0))**2 )
"""
ss.setIC(ic)

ss.updateVars()

ss.plot.figure(1)
ss.plot.clf()
ss.plot.contourf('phi',32)
ss.plot.figure(2)
ss.plot.clf()
ss.plot.contourf('ff',32)
plt.pause(.001)



            

