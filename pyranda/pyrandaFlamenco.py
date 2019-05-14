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
from .pyrandaPackage import pyrandaPackage
from . import parcop


class pyrandaFlamenco(pyrandaPackage):
    """
    Flamenco flux routines
    """
    def __init__(self,pysim):

        PackageName = 'Flamenco'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.order = 5   # Reconstruction scheme to use. 1=1st order,
        #                2=minmod, 3=van Leer, 4=superbee, 5=5th order MUSCL
        self.lowMach = 1  # Low mach correction 0=off, 1=on
        self.recon = 0    # Variable reconstruction (0=cons, 1=prim)
        self.eos = 1      # 1=ideal gas
        self.nSpec = 1    # Number of species

        self.grid = None  # Data for mesh
        
        
    def get_sMap(self):
        """
        String mappings for this package.  Packages added to the main
        pyranda object will check this map
        """
        sMap = {}
        sMap['flamenco.EulerFlux('] = "self.packages['Flamenco'].flux("
        self.sMap = sMap



    def flux(self,rhou,rhov,rhow,Et,rho,p,T,gamma):

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
        Grhou = self.pyranda.PyMPI.ghost( rhou, np=cHalo, perClip=False)
        Grhov = self.pyranda.PyMPI.ghost( rhov, np=cHalo, perClip=False)
        Grhow = self.pyranda.PyMPI.ghost( rhow, np=cHalo, perClip=False)
        Grho  = self.pyranda.PyMPI.ghost( rho , np=cHalo, perClip=False)
        GEt   = self.pyranda.PyMPI.ghost( Et,   np=cHalo, perClip=False)
        Gp    = self.pyranda.PyMPI.ghost( p,    np=cHalo, perClip=False)
        GT    = self.pyranda.PyMPI.ghost( T,    np=cHalo, perClip=False)
        Ggamma= self.pyranda.PyMPI.ghost( gamma,np=cHalo, perClip=False)
        
        
        # Get the mesh and get nodal data... (TODO)
        if self.grid == None:
            xc = self.pyranda.PyMPI.ghost(
                self.pyranda.mesh.coords[0].data,np=nHalo*2, perClip=False)
            yc = self.pyranda.PyMPI.ghost(
                self.pyranda.mesh.coords[1].data,np=nHalo*2, perClip=False)
            zc = self.pyranda.PyMPI.ghost(
                self.pyranda.mesh.coords[2].data,np=nHalo*2, perClip=False)

            # Fix boundary data for non-periodic meshes
            if self.pyranda.PyMPI.x1proc:
                xc[1,:,:] = 2.0*xc[2,:,:] - xc[3,:,:]
                xc[0,:,:] = 2.0*xc[1,:,:] - xc[2,:,:]
            if self.pyranda.PyMPI.xnproc:
                xc[-2,:,:] = 2.0*xc[-3,:,:] - xc[-4,:,:]
                xc[-1,:,:] = 2.0*xc[-2,:,:] - xc[-3,:,:]
            if self.pyranda.PyMPI.y1proc:
                yc[:,1,:] = 2.0*yc[:,2,:] - yc[:,3,:]
                yc[:,0,:] = 2.0*yc[:,1,:] - yc[:,2,:]
            if self.pyranda.PyMPI.ynproc:
                yc[:,-2,:] = 2.0*yc[:,-3,:] - yc[:,-4,:]
                yc[:,-1,:] = 2.0*yc[:,-2,:] - yc[:,-3,:]
            if self.pyranda.PyMPI.z1proc:
                zc[:,:,1] = 2.0*zc[:,:,2] - zc[:,:,3]
                zc[:,:,0] = 2.0*zc[:,:,1] - zc[:,:,2]
            if self.pyranda.PyMPI.znproc:
                zc[:,:,-2] = 2.0*zc[:,:,-3] - zc[:,:,-4]
                zc[:,:,-1] = 2.0*zc[:,:,-2] - zc[:,:,-3]

            # Get mean
            xn = numpy.zeros( (nx+1+nHalo*2,ny+1+nHalo*2,nz+1+nHalo*2))
            yn = numpy.zeros( (nx+1+nHalo*2,ny+1+nHalo*2,nz+1+nHalo*2))
            zn = numpy.zeros( (nx+1+nHalo*2,ny+1+nHalo*2,nz+1+nHalo*2))


            for i in range(0,xn.shape[0]):
                for j in range(0,yn.shape[1]):
                    for k in range(0,zn.shape[2]):
                        xn[i,j,k] = ( xc[i,j,k]       + xc[i+1,j,k] +
                                      xc[i,j+1,k]     + xc[i,j,k+1] +
                                      xc[i+1,j+1,k]   + xc[i+1,j,k+1] +
                                      xc[i+1,j+1,k+1] + xc[i,j+1,k+1] ) / 8.0
                        yn[i,j,k] = ( yc[i,j,k]       + yc[i+1,j,k] +
                                      yc[i,j+1,k]     + yc[i,j,k+1] +
                                      yc[i+1,j+1,k]   + yc[i+1,j,k+1] +
                                      yc[i+1,j+1,k+1] + yc[i,j+1,k+1] ) / 8.0
                        zn[i,j,k] = ( zc[i,j,k]       + zc[i+1,j,k] +
                                      zc[i,j+1,k]     + zc[i,j,k+1] +
                                      zc[i+1,j+1,k]   + zc[i+1,j,k+1] +
                                      zc[i+1,j+1,k+1] + zc[i,j+1,k+1] ) / 8.0

            self.grid = {}
            self.grid['xn'] = xn
            self.grid['yn'] = yn
            self.grid['zn'] = zn

        else:
            xn = self.grid['xn']
            yn = self.grid['yn']
            zn = self.grid['zn']
            
                        
        # Cv/Cp data arrays
        cv = numpy.zeros( (self.nSpec,2,5) )
        cp = numpy.zeros( (self.nSpec,2,5) )
        Rspec = numpy.zeros( self.nSpec )
        flux = parcop.parcop.flamencoflux(nx+1,ny+1,nz+1,xn,yn,zn,          # mesh
                                          nHalo,cHalo,                # ghosting
                                          self.order,self.recon,self.lowMach,
                                          Grhou,Grhov,Grhow,GEt,Grho,   # Cons
                                          Gp,Ggamma,GT)                 # con
                                          
                                
        
        Frhou = flux[0,cHalo:-cHalo,cHalo:-cHalo,cHalo:-cHalo]
        Frhov = flux[1,cHalo:-cHalo,cHalo:-cHalo,cHalo:-cHalo]
        Frhow = flux[2,cHalo:-cHalo,cHalo:-cHalo,cHalo:-cHalo]
        FEt   = flux[3,cHalo:-cHalo,cHalo:-cHalo,cHalo:-cHalo]
        Frho  = flux[4,cHalo:-cHalo,cHalo:-cHalo,cHalo:-cHalo]

        #import matplotlib.pyplot as plt
        #import pdb
        #pdb.set_trace()

        
        return [Frhou,Frhov,Frhow,FEt,Frho]
        

        
