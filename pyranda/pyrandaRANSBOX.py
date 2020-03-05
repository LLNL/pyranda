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

        # Inputs - auto from ransbox api
        inputs = ransbox.TurbSrcInputs()
        self.inputList = [i for i in dir(inputs) if "__" not in i]
        for ii in self.inputList:
            exec("self.%s = None" % ii)

        # Outputs - auto from ransbox api
        outputs = ransbox.TurbSrcTerms()
        self.outputList = [i for i in dir(outputs) if "__" not in i]
        for ii in self.outputList:
            exec("self.%s = None" % ii)
        
        
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

        # Setup inputs (if direct pointer, only do once?)
        inputs = ransbox.TurbSrcInputs()
        
        for ii in self.inputList:
            exec("val = self.%s" % ii )
            if val:
                exec("inputs.%s = self.toC( self.pysim.variables[ self.%s ].data )" % (ii,ii) )

            
        # Setup outputs
        outputs = ransbox.TurbSrcTerms()
        for ii in self.outputList:
            exec("val = self.%s" % ii )
            if val:
                exec("outputs.%s = self.toC( self.pysim.variables[ self.%s ].data )" % (ii,ii) )

        # Evaluate model
        a = RANSBox_EvaluateModel(self.model, inputs, outputs, 1, "cpu")


        for ii in self.outputList:
            exec("val = self.%s" % ii )
            if val:                    
                exec("self.pysim.variables[ self.%s ].data  = self.toF( outputs.%s )" % (ii,ii) )
                     
        
        
        
        
        

        
