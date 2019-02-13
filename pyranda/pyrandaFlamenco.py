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
from . import parcop


class pyrandaFlamenco(pyrandaPackage):
    """
    Flamenco flux routines
    """
    def __init__(self,name,pyranda):

        PackageName = 'Flamenco'
        pyrandaPackage.__init__(self,PackageName,pyranda)

        self.order = 5   # Reconstruction scheme to use. 1=1st order,
        #                2=minmod, 3=van Leer, 4=superbee, 5=5th order MUSCL
        self.lowMach = 1  # Low mach correction 0=off, 1=on
        self.recon = 0    # Variable reconstruction (0=cons, 1=prim)
        self.eos = 1      # 1=ideal gas
        self.nSpec = 1    # Number of species
        
    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        sMap['flamenco.EulerFlux('] = "self.packages['Flamenco'].flux("
        self.sMap = sMap



    def flux(self,rhou,rhov,rhow,Et,rho):

        # Grid size
        nx = self.pyranda.nx
        ny = self.pyranda.ny
        nz = self.pyranda.nz

        # Order sets ghost size
        cHalo = 3
        if self.order == 1:
            cHalo=1
        elif ( self.order >= 2 and self.order <=4 ):            
            cHalo=2
        elif (self.order == 5):
            cHalo=3
            
        # Number of halo nodes (always 1)
        nHalo = 1

        # Ghost the data and the mesh  TODO: may need BC aware ghosting
        gval = self.pyranda.PyMPI.ghost( val, np=3)
        # Get the mesh and get nodal data... (TODO)
        xc = self.pyranda.PyMPI.ghost(self.pyranda.mesh.coords[0].data,np=iNnHalo)
        yc = self.pyranda.PyMPI.ghost(self.pyranda.mesh.coords[1].data,np=iNnHalo)
        zc = self.pyranda.PyMPI.ghost(self.pyranda.mesh.coords[2].data,np=iNnHalo)


        # Cv/Cp data arrays
        cv = numpy.zeros( (self.nSpec,2,5) )
        cp = numpy.zeros( (self.nSpec,2,5) )
        Rspec = numpy.zeros( self.nSpec )
        flux = parcop.parcop.flamencoflux(nx,ny,nz,xc,yc,zc,          # mesh
                                          cHalo,nHalo,                # ghosting
                                          self.order,self.recon,self.lowMach,
                                          Grhou,Grhov,Grhow,GEt,Grho,   # Cons
                                          Gp,GT,Ggamma                   # con
                                          
                                
        

        return [Frhou,Frhov,Frhow,FEt,Frho]
        

        
