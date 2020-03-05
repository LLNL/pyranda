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

        
    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        sMap = ['RANS.eval()'] = "self.packages['RANSBOX'].evaluate()"
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
        inputs.den = self.toC( pysim.variables[ self.varMap['den'] ].data )
        inputs.tke = self.toC( pysim.variables[ self.varMap['tke'] ].data )
        inputs.Lt = self.toC( pysim.variables[ self.varMap['Lt'] ].data )
        inputs.nu = self.toC( pysim.variables[ self.varMap['nu'] ].data )

        # Setup outputs
        outputs = ransbox.TurbSrcTerms()
        outputs.vmut  = self.toC( pysim.variables[ self.varMap['vmut'] ].data )
        outputs.src_k = self.toC( pysim.variables[ self.varMap['src_k'] ].data )
        
        # Evaluate model
        a = RANSBox_EvaluateModel(self.model, inputs, outputs, 1, "cpu")

        # Reorder arrays (is this needed?)
        pysim.variables[ self.varMap['vmut'] ].data  = self.toF( outputs.vmut )
        pysim.variables[ self.varMap['src_k'] ].data = self.toF( outputs.src_k )
        
        
        
        
        

        
