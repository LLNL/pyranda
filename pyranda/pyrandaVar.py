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

class pyrandaVar:

    def __init__(self,name,kind,rank):

        self.name = name
        self.kind = kind
        self.rank = rank

        self.rhs = None   # RHS for algebraic closure
        self.dt  = None   # dt for temporal integration

        self.data = None  # Need to allocate
        self.PyMPI = None
        self.allocated = False
        
    def __allocate__(self,pympi):  # Make this private
        
        self.PyMPI = pympi

        # Inherit the mesh size
        if self.rank == 'scalar':
            self.data = pympi.emptyScalar()
        elif self.rank == 'vector':
            self.data = pympi.emptyVector()
        else:
            raise ValueError('Error: rank: %s not recognized' % self.rank)
        self.allocated = True


    def sum(self):
        
        if self.rank == 'scalar':
            return self.PyMPI.sum3D( self.data )
        elif self.rank == 'vector':
            return  [ self.PyMPI.sum3D( self.data[:,:,:,0] ) , 
                      self.PyMPI.sum3D( self.data[:,:,:,1] ) ,
                      self.PyMPI.sum3D( self.data[:,:,:,2] ) ]

    # Allow global access via indices
    def __getitem__(self, given):

        # Check length of indices
        slen = len( given )

        # Only 3D arguments
        if ( slen != 3 ):
            self.PyMPI.iprint("Error: 3 index ranges must be given")
            return None

        indices = [None]*3

        NX = [self.PyMPI.nx,self.PyMPI.ny,self.PyMPI.nz]
        for i in range(slen):
            
            if isinstance( given[i], slice):
                start = given[i].start
                stop = given[i].stop
                if start == None:
                    start = 0
                if stop == None:
                    stop = NX[i]
            else:
                start = given[i]
                stop = start + 1
            indices[i] = [ start, stop ]

        gdata = self.PyMPI.getIJK( self.data, indices[0],indices[1],indices[2] ) 
        return numpy.squeeze( gdata )

        

            
