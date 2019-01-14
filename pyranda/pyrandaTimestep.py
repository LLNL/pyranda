import numpy 
import re
import sys
import time
from .pyrandaPackage import pyrandaPackage


class pyrandaTimestep(pyrandaPackage):
    """
    Case physics package module for adding new physics packages to pyranda
    """
    def __init__(self,pysim):

        PackageName = 'Timestep'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.dx = pysim.mesh.d1
        self.dy = pysim.mesh.d2
        self.dz = pysim.mesh.d3
        self.pysim = pysim
        self.GridLen = pysim.mesh.GridLen
        

    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        sMap['dt.courant('] = "self.packages['Timestep'].courant("
        sMap['dt.diff('] = "self.packages['Timestep'].diff("
        self.sMap = sMap

    def courant(self,u,v,w,c):

        # Compute the dt for the courant limit
        if self.pysim.mesh.coordsys == 3:
            #dAdx = self.pysim.getVar("dAx")
            #dAdy = self.pysim.getVar("dAy")
            #dBdx = self.pysim.getVar("dBx")
            #dBdy = self.pysim.getVar("dBy")
            #vrate = ( numpy.abs(u*dAdx/self.pysim.dx + v*dAdy/self.pysim.dx) +
            #          numpy.abs(u*dBdx/self.pysim.dy + v*dBdy/self.pysim.dy) )
            umag = numpy.sqrt( u*u + v*v )
            vrate = umag / numpy.minimum( self.pysim.getVar('d1'),
                                          self.pysim.getVar('d2') )
                                      
            #vrate = ( numpy.abs(u) / self.dx + 
            #          numpy.abs(v) / self.dy + 
            #          numpy.abs(w) / self.dz )
        else:
            vrate = ( numpy.abs(u) / self.dx +
                      numpy.abs(v) / self.dy +
                      numpy.abs(w) / self.dz )
        
        
        crate = numpy.abs(c) / self.GridLen

        dt_max = 1.0 / self.pyranda.PyMPI.max3D(vrate + crate)

        return dt_max
        
        

    def diff(self,bulk,density):

        delta = self.GridLen
        drate = density * delta * delta / numpy.maximum( 1.0e-12, bulk )
        dt_max = self.pyranda.PyMPI.min3D( drate )

        return dt_max
        
        
