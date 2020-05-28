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
import ransbox
from ransbox import *
from .pyrandaPackage import pyrandaPackage

class pyrandaRANSBOX(pyrandaPackage):
    """
    RANSBOX package api
    """
    def __init__(self,pysim):

        PackageName = "RANSBOX"
        pyrandaPackage.__init__(self,PackageName,pysim)
        self.pysim = pysim                

        self.varMap = {}
        self.model  = None
        self.nx = pysim.PyMPI.ax
        self.ny = pysim.PyMPI.ay
        self.nz = pysim.PyMPI.az
        self.shape = (self.nx, self.ny, self.nz)
        
        # Inputs - auto from ransbox api
        inputs = ransbox.RANSBoxInputs()
        self.inputList = [i for i in dir(inputs) if "__" not in i]
        for ii in self.inputList:
            exec("self.%s = None" % ii)

        # Outputs - auto from ransbox api
        outputs = ransbox.RANSBoxOutputs()
        self.outputList = [i for i in dir(outputs) if "__" not in i]
        for ii in self.outputList:
            exec("self.%s = None" % ii)

        self.inputs  = inputs
        self.outputs = outputs

        self.allocated = False

    def allocate(self):

        # Allocate some space for RANSBox outputs
        npts = self.nx * self.ny * self.nz
        for ii in self.outputList:
            #exec("val = self.%s" % ii )
            #if val:
            exec("self.outputs.%s = numpy.zeros(npts)" % ii )

        self.allocated = True
                
        
    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        sMap['RANS.eval()'] = "self.packages['RANSBOX'].evaluate()"
        self.sMap = sMap

    # May require nothing for this stage
    def toC(self, val ):
        return val

    # Requires a re-ordering here.
    def toF(self, val ):
        # numpy
        val.shape = (self.nx, self.ny, self.nz)
        # val = val.reshape( self.nx, self.ny, self.nz, order = "F")  # Or maybe order = "C"
        return val
        
    def evaluate(self):

        if not self.allocated:
            self.allocate()
        
        # Setup inputs (if direct pointer, only do once?)
        #inputs = self.inputs
        
        for ii in self.inputList:
            exec("val = self.%s" % ii )
            if val:
                #exec("self.inputs.%s = self.toC( self.pysim.variables[ self.%s ].data )" % (ii,ii) )
                exec("self.inputs.%s = self.pysim.variables[ self.%s ].data " % (ii,ii) )

        # Evaluate model
        npts = self.nx * self.ny * self.nz
        RANSBox_EvaluateModel(self.model, self.inputs, self.outputs, npts, "cpu")

        #import pdb
        #pdb.set_trace()

        
        for ii in self.outputList:
            exec("val = self.%s" % ii )
            if val:                    
                #exec("self.pysim.variables[ self.%s ].data  = self.toF( self.outputs.%s )" % (ii,ii) )
                exec("self.pysim.variables[ self.%s ].data  = self.outputs.%s.reshape(self.shape)" % (ii,ii) )
                     
        
        
        
        
        

        
