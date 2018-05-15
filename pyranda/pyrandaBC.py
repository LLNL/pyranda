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
import re
import sys
import time
from pyrandaPackage import pyrandaPackage


class pyrandaBC(pyrandaPackage):

    def __init__(self,pysim):

        PackageName = 'BC'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.bcList = {}
        

    def get_sMap(self):
        sMap = {}
        sMap['bc.extrap('] =  "self.packages['BC'].extrap("
        sMap['bc.const(']  =  "self.packages['BC'].const("
        sMap['bc.exit(']   =  "self.packages['BC'].exitbc("
        self.sMap = sMap
        

    
        
    def extrap(self,var,direction):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.extrapolate( v , d )
            
        

    def extrapolate(self,var,direction):
        # Direction switch
        bcvar = None
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                self.pyranda.variables[var].data[0,:,:] = 2*self.pyranda.variables[var].data[1,:,:] - self.pyranda.variables[var].data[2,:,:] 
                #self.pyranda.variables[var].data[0,:,:] = self.pyranda.variables[var].data[1,:,:] 

        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                self.pyranda.variables[var].data[-1,:,:] = 2*self.pyranda.variables[var].data[-2,:,:] - self.pyranda.variables[var].data[-3,:,:]
                #self.pyranda.variables[var].data[-1,:,:] = self.pyranda.variables[var].data[-2,:,:]

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                self.pyranda.variables[var].data[:,0,:] = self.pyranda.variables[var].data[:,1,:]

        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                self.pyranda.variables[var].data[:,-1,:] = self.pyranda.variables[var].data[:,-2,:]

                
    def const(self,var,direction,val):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.constant( v , d , val)
                
                
    def constant(self,var,direction,val):
        
        # Direction switch
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                self.pyranda.variables[var].data[0,:,:] = val
                
            
        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                self.pyranda.variables[var].data[-1,:,:] = val
            

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                self.pyranda.variables[var].data[:,0,:] = val
                
            
        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                self.pyranda.variables[var].data[:,-1,:] = val
            
    def exitbc(self,var,direction,norm=False):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.exit_boundary( v , d , norm)

    def exit_boundary(self,var,direction,norm):


        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                u1 = self.pyranda.variables[var].data[0,:,:]
                u2 = self.pyranda.variables[var].data[1,:,:]
                u3 = self.pyranda.variables[var].data[2,:,:]                
                self.pyranda.variables[var].data[0,:,:] = self.BENO(u1,u2,u3,norm)
        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                u1 = self.pyranda.variables[var].data[-1,:,:]
                u2 = self.pyranda.variables[var].data[-2,:,:]
                u3 = self.pyranda.variables[var].data[-3,:,:]                
                self.pyranda.variables[var].data[-1,:,:] = self.BENO(u1,u2,u3,norm)



    def BENO(self,u1,u2,u3,norm):
        d1 =        u2     - u1
        d2 = (2.0*u2-u3) - u1

        beno = numpy.where( numpy.abs(d1) <= numpy.abs(d2), u1+d1, u1+d2)
        beno = numpy.where(d1*d2 <= 0.0, u1, beno)
              
        if norm:        
            beno = numpy.maximum( numpy.minimum( beno, numpy.maximum(0.0,u2) ), numpy.minimum( 0.0, u2) )

        return beno
