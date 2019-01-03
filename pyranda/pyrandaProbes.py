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
import sys
import time
from .pyrandaPackage import pyrandaPackage
from .pyrandaMisc import ipoint

class pyrandaProbes(pyrandaPackage):

    def __init__(self,pysim,x=None,y=None,z=None):

        PackageName = 'Probes'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.xpts = x
        self.ypts = y

        self.pysim = pysim
        # 3D not supported for now.
        self.zpts = None        
        
        self.values = None

        self.points = []

        xmesh = pysim.x().data
        ymesh = pysim.y().data
        for ix,iy in zip(self.xpts,self.ypts):
            self.points.append( ipoint(ix,iy,xmesh,ymesh) )
        
    def get(self,var,method='linear' ):

        # Get the pysim data
        values = []
        vals = self.pysim.variables[var].data[:,:,0]

        # Get interpolated values
        for pt in self.points:
            values.append( pt.interp( vals,method ) )

        # Communicate 
        Lvalues = numpy.array( values )
        Gvalues = self.pysim.PyMPI.comm.allreduce( Lvalues, op=MPI.SUM )
        
        self.values = Gvalues
        return self.values
        
    def get1(self,var,method='linear' ):

        # Get the pysim data
        values = []
        vals = self.pysim.variables[var].data[:,:,0]

        # Get interpolated values
        for pt in self.points:
            values.append( pt.interp1( vals,method ) )

        # Communicate 
        Lvalues = numpy.array( values )
        Gvalues = self.pysim.PyMPI.comm.allreduce( Lvalues, op=MPI.SUM )
        
        self.values = Gvalues
        return self.values
