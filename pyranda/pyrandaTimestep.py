import numpy 
import re
import sys
import time
from pyrandaPackage import pyrandaPackage


class pyrandaTimestep(pyrandaPackage):
    """
    Case physics package module for adding new physics packages to pyranda
    """
    def __init__(self,pysim):

        PackageName = 'Timestep'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.dx = None
        self.dy = None
        self.dz = None
        

    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        sMap['dt.courant('] = "self.packages['Timestep'].courant("
        sMap['dt.diff('] = "self.packages['Timestep'].diff("
        self.sMap = sMap

    def getDX(self):
        
        if not self.dx:
            x = self.pyranda.mesh.coords[0]
            y = self.pyranda.mesh.coords[1]
            z = self.pyranda.mesh.coords[2]

            try:
                self.dx = x[1,0,0] - x[0,0,0]
            except:
                self.dx = 1.0
            try:
                self.dy = y[0,1,0] - y[0,0,0]
            except:
                self.dy = 1.0
            try:
                self.dz = z[0,0,1] - z[0,0,0]
            except:
                self.dz = 1.0 

        return [self.dx,self.dy,self.dz]
                

    def courant(self,u,v,w,c):

        [dx,dy,dz] = self.getDX()

        # Compute the dt for the courant limit
        vrate = numpy.abs(u) / dx + numpy.abs(v) / dy + numpy.abs(w) / dz
        crate = numpy.abs(c) / dx

        dt_max = 1.0 / self.pyranda.PyMPI.max3D(vrate + crate)

        return dt_max
        
        

    def diff(self,bulk,density):

        [dx,dy,dz] = self.getDX()

        delta = (dx*dy*dz)**(1./3.)
        drate = density * delta * delta / numpy.maximum( 1.0e-12, bulk )

        dt_max = self.pyranda.PyMPI.min3D( drate )

        return dt_max
        
        
