from __future__ import print_function
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams

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


# Define the equations of motion
eom = """
# Heat equation
ddt(:T:)  =  lap(:T:) + pi*(2.0 - lap(:psi:))*cos(pi*:psi:) + pi**2*(:psiX:**2 + :psiY:**2)*sin(pi*:psi:)
:error: = abs( :T: - :bb: ) * where( :phi: < 0.0, 0.0 , 1.0)
:psi: = 9.0*( (meshx-pi)/4.0 + 1)**2 + 4.0*( (meshy-pi)/4.0 + 1)**2 + 2.*simtime
[:psiX:, :psiY:, :psiZ: ] = grad( :psi: )
:bb: = sin( pi * :psi: )
:T: = ibmC( :T: , :phi:, [:gx:,:gy:,:gz:] , :bb:)
"""


# Initialize variables
ic = """
xnn = meshx[-1,0,0]
:psi: = 9.0*( (meshx-pi)/4.0 + 1)**2 + 4.0*( (meshy-pi)/4.0 + 1)**2 + 2.*simtime
:T: = sin( pi * :psi: )
:bb: = :T: * 1.0
rad =  sqrt( (meshx-pi)**2 + (meshy-pi)**2 )
theta = numpy.arctan2(meshy-pi,meshx-pi)
xth = (10.*sin(2*theta)**2 + 3.*cos(2.*theta)**3 + 40 )* cos(theta)/20. + pi
yth = (10.*sin(2*theta)**2 + 3.*cos(2.*theta)**3 + 40 )* sin(theta)/20. + pi
rth = sqrt( (xth-pi)**2 + (yth-pi)**2 )
:phi:  = rth - rad 
[:gx:,:gy:,:gz:] = grad( :phi: )
"""

rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 20
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True



# Initialize a simulation object on a mesh
def solve( Npts,test=True):

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
    
    ss = pyrandaSim('heat_equation',mesh_options,silent=test)
    ss.addPackage( pyrandaBC(ss) )
    ss.addPackage( pyrandaIBM(ss) )
    ss.addPackage( pyrandaTimestep(ss) )
    ss.EOM(eom)
    ss.setIC(ic)
    x  = ss.mesh.coords[0].data
    xx =  ss.PyMPI.zbar( x )

    # Time step size
    dt_max = (L / ss.mesh.nn[0])**2 * .10
    tt = 1.0 #dt_max * 2000 


    # Main time loop for physics
    dt = dt_max
    cnt = 1
    time = 0.0
    viz = True
    while tt > time:

        time = ss.rk4(time,dt)
        dt = min(dt_max, (tt - time) )

        if not test:
            ss.iprint("%s -- %s" % (cnt,time)  )

        cnt += 1
        if viz and (not test):

            # Plot animation of advection
            if (ss.PyMPI.master and (cnt%200 == 0)):
                ss.plot.figure(1)
                ss.plot.clf()

                #ss.plot.contourf('T',64,cmap=cm.jet)
                ss.plot.contourf('error',64,cmap=cm.jet)
                ss.plot.contour('phi',[0.0])

                ss.plot.figure(2)
                ss.plot.clf()
                ss.plot.contourf('bb',64,cmap=cm.jet)
                ss.plot.contour('phi',[0.0])
                plt.pause(.001)

                plt.xlabel(r"$x$") 
                plt.ylabel(r"$y$")
                #import pdb
                #pdb.set_trace()

    # Plot animation of advection
    error = ss.PyMPI.zbar( ss.variables['error'].data )

    L2     = numpy.sum( error  ) / Npts**2
    Linfty = numpy.max( error  )

    print( "L_infinity = %s" % Linfty )
    print( "L2        = %s" % L2 ) 

    return [ss,Linfty,L2]


# Results for quadratic BCs (constant and linear are exact)

RES = [50,75,100,150,200]
L2 = []
Li = []
dx = []
for res in RES:
    [ss,myLi,myL2] = solve(res,test=True)
    L2.append(myL2)
    Li.append(myLi)
    dx.append(2.0*numpy.pi / float(res))


LW = 2.5

# 8th order fiducial
DX = numpy.array(dx)
plt.loglog(1.0/DX,Li,'bo', linewidth=LW)
plt.loglog(1.0/DX,DX**4*1.0e1,'k--')
plt.loglog(1.0/DX,DX**10 * 1.0e8,'k--')
#plt.loglog(1.0/DX,DX**2*1.0e0,'k--')

plt.tight_layout(pad=3.0)
plt.xlabel(r"$\frac{1}{\Delta x}$") 
plt.ylabel(r"$L_\infty$ Error")

floc = "/Users/olson45/Documents/Conference Travel/ISSW32_2019/ISSW32_bjo/figures"
import os

plt.savefig(os.path.join(floc,'heat_2D_conv.pdf'))
plt.show()


#python ISSW_heatSIN.py 32 1
#L1 = 3.08891435063
#L2 = 0.109368403694

#python ISSW_heatSIN.py 64 1
#L1 = 0.013809734169
#L2 = 0.000330344858814

#python ISSW_heatSIN.py 128 1
#L1 = 1.9008337412e-05
#L2 = 4.37165809992e-06

#python ISSW_heatSIN.py 256 1
#L1 = 1.70971144206e-05
#L2 = 4.41830383029e-06
