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
import scipy.interpolate as SI
import scipy
import re
import sys
import time
from .pyrandaPackage import pyrandaPackage

class pyrandaProbes(pyrandaPackage):

    def __init__(self,pysim,x=None,y=None,z=None):

        PackageName = 'IBM'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.xpts = x
        self.ypts = y

        # 3D not supported for now.
        self.zpts = None        
        
        self.values = None

        self.points = []

        xmesh = pymsim.grid('x').data
        ymesh = pymsim.grid('y').data
        for ix,iy in zip(self.xpts,self.ypts):
            self.points.append( ipoint(ix,iy,xmesh,ymesh) )
        
    def get(self,var ):

        self.values = []
        vals = pysim.variables[var].data
        
        for pt in self.points:
            self.values.append( pt.interp( vals ) )
        
        
            
        
class ipoint:
    """
    Interpolated object to rapidly get values
      Note: grids should be ghosted in the +i,j+ directions
    """
    def __init__(self,xpt,ypt,xgrid,ygrid):
        

        self.value = 0.0
        self.onProc = False
        
        self.x = xpt
        self.y = ypt

        # Bounding box for determining proc location
        onX = ( (self.x > xgrid.min() ) and  (self.x < xgrid.max() ) )
        onY = ( (self.y > ygrid.min() ) and  (self.y < ygrid.max() ) )        
        if onX and onY:
            self.onProc = True
        else:
            return
            
        self.ii1 = -1
        self.iin = -1

        self.jj1 = -1
        self.jjn = -1

        [self.nx,self.ny] = xgrid.shape
               
        self.weight  = []

        x_y = numpy.dstack([xgrid.ravel(),ygrid.ravel()])[0]

        pts = numpy.zeros( (2,1) )
        pts[0] = self.x*1.0
        pts[1] = self.y*1.0
        points = list(pts.transpose())

        mytree = scipy.spatial.cKDTree(x_y)
        dist, indexes = mytree.query(points)

        ii = int(indexes[0]/self.ny)
        jj = indexes[0]%self.ny

        x0 = numpy.abs(xgrid[ii,jj] - self.x)
        y0 = numpy.abs(ygrid[ii,jj] - self.y)

        self.ii0 = ii
        self.jj0 = jj

        self.xlocal = numpy.array([ [xgrid[ii-1,jj-1],xgrid[ii,jj-1],xgrid[ii+1,jj-1]],
                                    [xgrid[ii-1,jj  ],xgrid[ii,jj  ],xgrid[ii+1,jj  ]],
                                    [xgrid[ii-1,jj+1],xgrid[ii,jj+1],xgrid[ii+1,jj+1]] ])

        self.ylocal = numpy.array([ [ygrid[ii-1,jj-1],ygrid[ii,jj-1],ygrid[ii+1,jj-1]],
                                    [ygrid[ii-1,jj  ],ygrid[ii,jj  ],ygrid[ii+1,jj  ]],
                                    [ygrid[ii-1,jj+1],ygrid[ii,jj+1],ygrid[ii+1,jj+1]] ])

        

    def interp(self,vals,imethod='linear'):


        if (not self.onProc):
            return 0.0
        
        ii = self.ii0
        jj = self.jj0
        vlocal = numpy.array([ [vals[ii-1,jj-1],vals[ii,jj-1],vals[ii+1,jj-1]],
                             [vals[ii-1,jj  ],vals[ii,jj  ],vals[ii+1,jj  ]],
                             [vals[ii-1,jj+1],vals[ii,jj+1],vals[ii+1,jj+1]] ])
        
        pts = [ numpy.array([self.x]), numpy.array([self.y]) ]
        dc = self.getPoints(self.xlocal,self.ylocal,vlocal,pts,imethod=imethod)[0]

        return dc


    def getPoints(x,y,f,pts,imethod='linear'):
        xx = x.ravel()
        yy = y.ravel()
        ff = f.ravel()
        xpts = pts[0].ravel()
        ypts = pts[1].ravel()

        dc = SI.griddata(npy.array([xx,yy]).T,ff,npy.array([xpts,ypts]).T,method=imethod)
    
        return dc



"""
from mpi4py import MPI
import numpy
import scipy.interpolate as SI
import scipy
import re
import sys
import time


pts = numpy.zeros( (2,1) )

pts[0] = 15
pts[1] = 15

points = list(pts.transpose())

[xgrid,ygrid] = numpy.meshgrid( numpy.linspace(0,10,100), numpy.linspace(0,10,100) )

x_y = numpy.dstack([xgrid.ravel(),ygrid.ravel()])[0]
mytree = scipy.spatial.cKDTree(x_y)
dist, indexes = mytree.query(points)

ii = int(indexes[0]/100)
jj = indexes[0]%100


"""
