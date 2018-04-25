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
class pyrandaVar:

    def __init__(self,name,kind,rank):

        self.name = name
        self.kind = kind
        self.rank = rank

        self.rhs = None   # RHS for algebraic closure
        self.dt  = None   # dt for temporal integration

        self.data = None  # Need to allocate
        self.PyMPI = None

    def __allocate__(self,pympi):  # Make this private
        
        self.PyMPI = pympi

        # Inherit the mesh size
        if self.rank == 'scalar':
            self.data = pympi.emptyScalar()
        elif self.rank == 'vector':
            self.data = pympi.emptyVector()
        else:
            raise ValueError('Error: rank: %s not recognized' % self.rank)



    def sum(self):
        
        if self.rank == 'scalar':
            return self.PyMPI.sum3D( self.data )
        elif self.rank == 'vector':
            return  [ self.PyMPI.sum3D( self.data[:,:,:,0] ) , 
                      self.PyMPI.sum3D( self.data[:,:,:,1] ) ,
                      self.PyMPI.sum3D( self.data[:,:,:,2] ) ]
