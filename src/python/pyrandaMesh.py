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
import numpy 
from pyrandaUtils import *


class pyrandaMesh:

    def __init__(self):#,mesh_options):

        self.name = 'base'
        self.kind = None #mesh_options['type']
        self.options = {} #None #mesh_options
        self.dims = 0 #mesh_options['dim']
        self.PyMPI = None
        self.x = None
        self.y = None
        self.z = None
        self.coords = None
        self.shape = None
        #self.makeMesh()

    def makeMeshStr(self,str_mesh):

        mesh_options = {}
        meshStrLines = splitLines( str_mesh ) # Split into lines

        
        self.options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
        self.options['xn'] = [ 1   , 1    ,  1   ]
        self.options['nn'] = [ 1   , 1    ,  1   ]
        self.options['type'] = 'cartesian'
        self.options['periodic'] = numpy.array([False, False, False])
        self.options['dim'] = 1
        
        self.get_sMap()

        for msl in meshStrLines:
            exec( fortran3d(msl,self.sMap) )
        
        
    def set_options(self,ind,x1,xn,nn,periodic=False):

        self.options['x1'][ind] = x1
        self.options['xn'][ind] = xn 
        self.options['nn'][ind] = nn
        self.options['periodic'][ind] = periodic 
        
    def get_sMap(self):
        sMap = {}
        sMap['xdom=('] = "self.set_options(0,"
        sMap['ydom=('] = "self.set_options(1,"
        sMap['zdom=('] = "self.set_options(2,"
        self.sMap = sMap

    def makeMesh(self):


        options = self.options
        if self.kind == 'cartesian':

            x1 = options['x1']
            xn = options['xn']
            nn = options['nn']

            self.nn = nn    

            chunk_lo = self.PyMPI.chunk_3d_lo
            chunk_hi = self.PyMPI.chunk_3d_hi
            x = numpy.linspace(x1[0],xn[0],num=nn[0]+1)[chunk_lo[0]:chunk_hi[0]+1]
            y = numpy.linspace(x1[1],xn[1],num=nn[1]+1)[chunk_lo[1]:chunk_hi[1]+1]
            z = numpy.linspace(x1[2],xn[2],num=nn[2]+1)[chunk_lo[2]:chunk_hi[2]+1]

            x = numpy.linspace(x1[0],xn[0],num=nn[0])[chunk_lo[0]:chunk_hi[0]+1]
            y = numpy.linspace(x1[1],xn[1],num=nn[1])[chunk_lo[1]:chunk_hi[1]+1]
            z = numpy.linspace(x1[2],xn[2],num=nn[2])[chunk_lo[2]:chunk_hi[2]+1]
            
            x, y, z = numpy.meshgrid(x, y, z, indexing='ij')

            self.coords = [0]*3
            self.coords[0] = numpy.asfortranarray( x )
            self.coords[1] = numpy.asfortranarray( y )
            self.coords[2] = numpy.asfortranarray( z )

            self.shape = list(self.coords[0].shape)

