from __future__ import print_function
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM

# Try to get args
try:
    Npts = int(sys.argv[1])
except:
    Npts = 100

try:
    test = bool(sys.argv[2])
except:
    test = False

## Define a mesh
L = numpy.pi * 2.0  
Lp = L * (Npts-1.0) / Npts


mesh_options = {}
mesh_options['coordsys'] = 0
mesh_options['periodic'] = numpy.array([False, True, True])
mesh_options['dim'] = 1
mesh_options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
mesh_options['xn'] = [ Lp  , Lp   ,  Lp  ]
mesh_options['nn'] = [ Npts, Npts ,  1   ]




# Define the equations of motion
eom = """
# Heat equation
ddt(:T:)  =  :c: * lap(:T:)
#:T: = ibmS( :T: , :phi:, [:gx:,:gy:,:gz:] )
:T: = ibmC( :T: , :phi:, [:gx:,:gy:,:gz:] , :bb:)
:error: = abs( :T: - :bb: ) * where( :phi: < 0.0, 0.0 , 1.0)
"""


# Initialize variables
ic = """
xnn = meshx[-1,0,0]
:T: = 1.0 + 1.0*(xnn - meshx)**func/xnn**func
:bb: = :T: * 1.0
:c:   = 1.0
rad =  sqrt( (meshx-pi)**2 + (meshy-pi)**2 )
:phi:  = 1.5 - rad 
[:gx:,:gy:,:gz:] = grad( :phi: )
"""





# Initialize a simulation object on a mesh
def solve( Npts, func=0,viz=False,test=True):

    ss = pyrandaSim('heat_equation',mesh_options)
    ss.addPackage( pyrandaBC(ss) )
    ss.addPackage( pyrandaIBM(ss) )
    ss.addPackage( pyrandaTimestep(ss) )
    ss.EOM(eom)
    ss.setIC(ic.replace('func',str(func)))
    x  = ss.mesh.coords[0].data
    xx =  ss.PyMPI.zbar( x )

    # Time step size
    dt_max = L / ss.mesh.nn[0] * .005
    tt = dt_max * 2000 


    # Main time loop for physics
    dt = dt_max
    cnt = 1
    time = 0.0
    #viz = True
    while tt > time:

        time = ss.rk4(time,dt)
        dt = min(dt_max, (tt - time) )

        if not test:
            ss.iprint("%s -- %s" % (cnt,time)  )
            
        cnt += 1
        if viz:

            # Plot animation of advection
            if (ss.PyMPI.master and (cnt%50 == 0)):
                plt.figure(1)
                plt.clf()

                #ss.plot.contourf('T',64,cmap=cm.jet)
                ss.plot.contourf('error',64,cmap=cm.jet)
                ss.plot.contour('phi',[0.0])
                plt.pause(.001)

                import pdb
                pdb.set_trace()

    # Plot animation of advection
    error = ss.PyMPI.zbar( ss.variables['error'].data )
    
    L2 = numpy.sum( error  ) / Npts**2
    L1 = numpy.max( error  )
    
    print( "L1 = %s" % L1 )
    print( "L2 = %s" % L2 ) 

    return ss


# Results for quadratic BCs (constant and linear are exact)
#solve(50,func=2)
#0.160014424928
#solve(100,func=2)
#0.0400036062319
#solve(200,func=2)
#0.010000901558
#solve(400,func=2)
#0.0025002253895

# Second order
