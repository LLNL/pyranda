################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
################################################################################
from mpi4py import MPI
import numpy 
import re
import sys,os
import time
import matplotlib.pyplot as plt
from matplotlib import cm

from .pyrandaMPI   import pyrandaMPI
from .pyrandaVar   import pyrandaVar
from .pyrandaEq    import pyrandaEq
from .pyrandaMesh  import pyrandaMesh
from .pyrandaIO    import pyrandaIO
from .pyrandaUtils import *
from .pyrandaTex   import pyrandaTex



class pyrandaPlot:

    def __init__(self,pyranda):

        self.pyranda = pyranda  # point to main obj

        
    # Some master overloaded functions
    def figure(self,fignum):
        if self.pyranda.PyMPI.master:
            plt.figure(fignum)

    def show(self):
        if self.pyranda.PyMPI.master:
            plt.show()

    def clf(self):
        if self.pyranda.PyMPI.master:
            plt.clf()

    def title(self,name):
        if self.pyranda.PyMPI.master:
            plt.title(name)

    def pause(self):
        if self.pyranda.PyMPI.master:
            plt.pause()

    def contour(self,var,levels, slice3d='k=0',**kwargs):
        self.contourf(var , levels, slice3d=slice3d,filled=False,**kwargs)
            
    def contourf(self, var , levels, slice3d='k=0',filled=True,**kwargs):
        """
        Plots 2d contour plot.  Only works for 2d and 3d problems.
        2d - figure our directions and plot
        3d - assumes slice at k=0 unless specified
        """
        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz
        # 1D - Reject
        if (ny == nz == 1) or (nx == nz == 1) or (ny == nx == 1):
            self.pyranda.iprint("pyranda contourf only works for 2D & 3D problems")

        x = self.pyranda.mesh.coords[0]
        y = self.pyranda.mesh.coords[1]
        vv = self.pyranda.variables[var].data

        
        # 2D - figure out directions
        if (nz == 1):
            xdata = self.pyranda.PyMPI.zbar( x )
            ydata = self.pyranda.PyMPI.zbar( y )
            vdata = self.pyranda.PyMPI.zbar( vv )


        # 3D - pick slice
        if not (nx==1 or ny==1 or nz==1):

            # get directions - k slice
            if 'k' in slice3d:
                exec(slice3d)
                xdata = self.pyranda.PyMPI.getIJK(x,[0,nx],[0,ny],[k,k+1])[:,:,0]
                ydata = self.pyranda.PyMPI.getIJK(y,[0,nx],[0,ny],[k,k+1])[:,:,0]
                vdata = self.pyranda.PyMPI.getIJK(vv,[0,nx],[0,ny],[k,k+1])[:,:,0]
                
        if self.pyranda.PyMPI.master:
            if filled:
                plt.contourf(xdata,ydata,vdata,levels,**kwargs)
                plt.colorbar()
            else:
                if not (type(levels) == type([]) ):
                    levels = numpy.linspace(vdata.min(),vdata.max(),levels)
                plt.contour(xdata,ydata,vdata,levels,**kwargs)
            plt.axis('equal')
            plt.pause(.01)

                

    def showGrid(self):

        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz
        
        # 2D - figure out directions
        if (nz == 1):
            x  = self.pyranda.mesh.coords[0]
            y  = self.pyranda.mesh.coords[1]
            xx = self.pyranda.PyMPI.zbar( x )
            yy = self.pyranda.PyMPI.zbar( y )
        
        if self.pyranda.PyMPI.master:
            plt.plot(xx, yy, 'k-', lw=0.5, alpha=0.5)
            plt.plot(xx.T, yy.T, 'k-', lw=0.5, alpha=0.5)
            plt.axis('equal')
            plt.pause(.01)
