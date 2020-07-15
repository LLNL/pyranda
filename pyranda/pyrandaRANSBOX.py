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
try:
    import ransbox
    from ransbox import *
except:
    pass
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

        # Available Models
        self.modPrefix = "RANSBox_ModelType_"
        self.availableModels = self.getModelsList()

        # Available Coordsys
        self.coordPrefix = "RANSBox_Coord_"
        self.availableCoords = self.getCoordsList()

        
        self.inputs  = inputs
        self.outputs = outputs
        self.allocated = False

    def createModel(self,modelName,coordName):
        """Given a model name, create the model"""

        # Check that model is valid
        modelName = modelName.replace( self.modPrefix, "")
        RBmodelName = self.modPrefix + modelName
        valid = False
        for mt in self.availableModels:
            if RBmodelName in mt:
                valid = True

        if not valid:
            self.pysim.iprint("Selected model: %s -- (%s) is Not valid" % (modelName,RBmodelName))
            self.pysim.iprint("Valid models are (select either string):")
            self.printModels()
            return

        model = eval(RBmodelName)

        # Check that coords is valid
        coordName = coordName.replace( self.coordPrefix, "")
        RBcoordName = self.coordPrefix + coordName
        valid = False
        for mt in self.availableCoords:
            if RBcoordName in mt:
                valid = True

        if not valid:
            self.pysim.iprint("Selected coord system: %s --- (%s) is Not valid" % (coordName,RBcoordName))
            self.pysim.iprint("Valid coordinate are (select either string):")
            self.printCoords()
            return

        coord = eval(RBcoordName)
        
        # Make the model
        self.model = RANSBox_CreateModel(model,coord)
        
        

    def getCoordsList(self):
        """ Get available RANSBOX coordinate types"""
        availableCoords = []
        for d in dir(ransbox):
            if self.coordPrefix in d:
                availableCoords.append(d)
        return availableCoords

    def printCoords(self):
        """ Prints available RANSBOX coords"""
        for d in self.availableCoords:
            self.pysim.iprint("%s --- %s" % (d.replace(self.coordPrefix,''),d))
                                                        
    def getModelsList(self):
        """ Get available RANSBOX models"""
        availableModels = []
        for d in dir(ransbox):
            if self.modPrefix in d:
                availableModels.append(d)
        return availableModels

    def printModels(self):
        """ Prints available RANSBOX models"""
        for d in self.availableModels:
            self.pysim.iprint("%s --- %s" % (d.replace(self.modPrefix,''),d))

    def getSelfSimilarity(self):
        """Get an object for setting SS parameters"""
        parms = ransbox.RANSBoxSimilarityParams()
        err = ransbox.RANSBox_GetSelfSimilarityParams(self.model, parms)
        if ( err == 1 ):
            self.pysim.iprint("Error getting RANSBOX self similarity object")
            raise ValueError()
        return parms

    def setSelfSimilarity(self,parms):
        """Set the parms object as SS parameters"""
        err = ransbox.RANSBox_SetSelfSimilarModelCoeffs(self.model,parms)
        if ( err == 1 ):
            self.pysim.iprint("Error setting RANSBOX self similarity object")
            raise ValueError()

    def setModelCoeff(self,name,val):
        err = ransbox.RANSBox_SetModelCoeffByName(self.model,name,val)
        if ( err == 1 ):
            self.pysim.iprint("Error setting RANSBOX Coefficient: %s not valid" % name)
            raise ValueError()
        
        
                
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
                     
        
        
        
        
        

        
