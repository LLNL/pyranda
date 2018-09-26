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

        self.dx = pysim.dx
        self.dy = pysim.dy
        self.dz = pysim.dz

        self.GridLen = pysim.PyMPI.GridLen
        

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
        vrate = numpy.abs(u) / self.dx + numpy.abs(v) / self.dy + numpy.abs(w) / self.dz
        crate = numpy.abs(c) / self.GridLen

        dt_max = 1.0 / self.pyranda.PyMPI.max3D(vrate + crate)

        return dt_max
        
        

    def diff(self,bulk,density):

        delta = self.GridLen
        drate = density * delta * delta / numpy.maximum( 1.0e-12, bulk )
        dt_max = self.pyranda.PyMPI.min3D( drate )

        return dt_max
        
        
