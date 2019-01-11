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
from mpi4py import MPI
import numpy 
import re
import sys
import time
from .pyrandaPackage import pyrandaPackage


class pyrandaBC(pyrandaPackage):

    def __init__(self,pysim):

        PackageName = 'TBL'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.bcList = {}

        self.BCdata = {}


        # Persistent 2D arrays (mesh size)
        self.A11 = None
        self.A22 = None
        self.A33 = None
        self.A12 = None
        self.Umean = None
        self.Tmean = None
        self.RHO_u  = None
        self.RHO_v = None
        self.RHO_w = None
        self.buI = None
        self.bvI = None
        self.bwI = None
        self.buO = None
        self.bvO = None
        self.bwO = None
        self.y_r = None

        # Global array
        self.y_G = None

        # Integers
        self.UIx = 0
        self.UIz = 0
        self.VIx = 0
        self.VIz = 0
        self.WIx = 0
        self.WIz = 0

        self.UIy = [0,0]
        self.VIy = [0,0]
        self.WIy = [0,0]

        self.Nbuff = 0


        

    def get_sMap(self):
        sMap = {}
        sMap['TBL.inflow('] =  "self.packages['TBL'].inflow("
        sMap['TBL.setup('] =  "self.packages['TBL'].setup("
        self.sMap = sMap
        

    def setup(self,TBLfile,flowdir,walldir):

        # Flow in -x- wall in -y-
        aF = self.PyMPI.ax
        aW = self.PyMPI.ay
        aT = self.PyMPI.az

        nW = self.PyMPI.ny
        
        self.Umean = numpy.zeros( (aW,aT) )
        self.Tmean = numpy.zeros( (aW,aT) )
        self.A11 = numpy.zeros( (aW,aT) )
        self.A22 = numpy.zeros( (aW,aT) )
        self.A33 = numpy.zeros( (aW,aT) )
        self.A12 = numpy.zeros( (aW,aT) )

        self.RHO_u = numpy.zeros( (aW,aT) )
        self.RHO_v = numpy.zeros( (aW,aT) )
        self.RHO_w = numpy.zeros( (aW,aT) )

        self.y_G = numpy.zeros( (nW) )
        self.y_r = numpy.zeros( (aW) )


        tmp = self.pyranda.PyMPI.emptyScalar()
        if self.pyranda.PyMPI.x1proc:  # Flow dir
            tmp = self.pyranda.mesh.coords.data[1]  # Wall dir
            tmp /= tmp.size[0]
            tmp /= tmp.size[2]
        tmp1 = numpy.sum( tmp, (0,2) )
            
        self.y_G = numpy.concatenate( self.pyranda.PyMPI.ycom.allgather( tmp1 ) )
        self.y_r = self.y_G[0:aW]
            

        # Read in BL data from file
        #   and place on 2D planes Umean/Tmean/UUmean/VVmean/WWmean/UVmean
        self.A11 = numpy.sqrt( UUmean )
        self.A12 = UVmean / A11
        numpy.where( self.A11 == 0.0, self.A12 = 0.0 )
        self.A22 = numpy.sqrt( VVmean - self.A12**2 )
        self.A33 = numpy.sqrt( WWmean )
        
        # Calculate the Displacement thickness
        del_star = 0.0
        for i in range(ny/2):
            del_star = del_star + (1.0 - Uprof[i])* (y_G[i+1]-y_G[i])
            del_star = del_star * del_BL  # The U prof was non-dim by del_BL

        # Non-dimensionalize the local wall normal length scale for filter size selection
        y_r = y_r / del_star

        # Setup the 2d filters ( 6 total.  Inner (u,v,w) and outer (u,v,w) )
        self.UIx = 10
        self.UIy = [ 20, 35 ]
        self.UIz = 20
        
        self.VIx = 4
        self.VIy = [ 25, 45 ]
        self.VIz = 20
      
        self.WIx = 4
        self.WIy = [ 15, 20 ]
        self.WIz = 30

        Nbuff = 45 * 2

        self.simtime_old = self.pyranda.simtime
        self.tauX = float(UIx)*del_star / U_in

        # Inner coefficients
        self.buI = setFiltCoeffs(UIy[0],UIz)
        self.bvI = setFiltCoeffs(VIy[0],VIz)
        self.bwI = setFiltCoeffs(WIy[0],WIz)

        # Outer coefficients
        self.buO = setFiltCoeffs(UIy[1],UIz)
        self.bvO = setFiltCoeffs(VIy[1],VIz)
        self.bwO = setFiltCoeffs(WIy[1],WIz)

        # Set up a restart directory
        # // TODO //

        
        # Initialize the temporal averages here.  They will be 
        #   over written on restart read. (done)

        # TODO: READ IN OLD DF FILES TO PICK UP THE RHO_u, RHO_v. RHO_w


        # DONE SETUP Digitial filtering
        

    def DFinflow(self):

        # Call fortran each time step

        
        rands = self.get_rands( ny, nz, Nbuff )
        # make parcop operator for this...
        rands = parcop.get_rands_normal( ny, nz, Nbuff ) 


        vU = parcop.filtRands(self.UIz,self.UIy[0],self.UIy[1],
                              self.Nbuff,self.buI,self.buO,rands[0,:,:])
        vV = parcop.filtRands(self.VIz,self.VIy[0],self.VIy[1],
                              self.Nbuff,self.bvI,self.bvO,rands[1,:,:])
        vW = parcop.filtRands(self.WIz,self.WIy[0],self.WIy[1],
                              self.Nbuff,self.bwI,self.bwO,rands[2,:,:])
                         

        # Restart IO


        # Get the updated rho_k
        # Time avergaging coefficients and quantities/fluctuations
        # Need to write a restart file with averages stored
        t_nm1 = self.simtime_old*1.0
        t_n = self.pyranda.simtime
        self.simtime_old = t_n
        dt_n = t_n - t_nm1
        EXPt = numpy.exp( - numpy.pi * dt_n / self.tauX )
        self.RHO_u += numpy.sqrt( EXPt ) + vU * numpy.sqrt( 1.0- EXPt )
        self.RHO_v += numpy.sqrt( EXPt ) + vV * numpy.sqrt( 1.0- EXPt )
        self.RHO_w += numpy.sqrt( EXPt ) + vW * numpy.sqrt( 1.0- EXPt )

        # Add the perturbations to mean with given 2 point correlations
        uinlet = self.Umean +  self.A11 * self.RHO_u
        vinlet =               self.A12 * self.RHO_u + self.A22 * self.RHO_v  
        winlet =               self.A33 * self.RHO_w

        # Get the temperature perturbation using strong Reynolds Analogy (SRA)
        MaSq = Mach**2
        Tinlet = self.Tmean + self.Tmean * ( -gm1*MaSq * (uinlet - self.Umean) / U_in )
        
        # Add to inlet
        if self.pyranda.PyMPI.x1proc:
            self.pyranda.variables[self.u].data[0,:,:] = uinlet
            self.pyranda.variables[self.v].data[0,:,:] = vinlet
            self.pyranda.variables[self.w].data[0,:,:] = winlet
            self.pyranda.variables[self.T].data[0,:,:] = Tinlet

        # Reconcile EOS outside of this loop/BC call
            

        
    def inflow(self,var,direction,field):

        val = field
        
        # Direction switch
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                self.pyranda.variables[var].data[0,:,:] = val[0,:,:]
                
            
        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                self.pyranda.variables[var].data[-1,:,:] = val[-1,:,:]
            

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                self.pyranda.variables[var].data[:,0,:] = val[:,0,:]
                
            
        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                self.pyranda.variables[var].data[:,-1,:] = val[:,-1,:]

        

                 

    def get_rands(self,ny,nz,Nbuff):

        # Get a common seed
        seed = 0
        if self.pyranda.PyPMI.master:
            seed = int(numpy.random.rand()*2**32)
            seed = self.pyranda.PyMPI.xyzcom.allreduce(seed, op=MPI.SUM)

        numpy.random.seed( seed )
        
        return numpy.random.rand(4,ny+Nbuff,nz+Nbuff)
        
        

                
def setFiltCoeffs(self,N1,N2):


    bij = numpy.zeros( (N1,N2) )
    
    den1 = calcDen(N1)
    den2 = calcDen(N2)

    N1d2 = float(N1)/2.0
    N2d2 = float(N2)/2.0

    for i in range(1,2*N1+2):
        for j in range(1,2*N2+2):
            k1 = abs(i-(N1+1))
            k2 = abs(j-(N2+1))
            bij[i,j] = numpy.exp(-numpy.pi*float(k1)/N1d2)/den1 * numpy.exp(-pi*float(k2)/N2d2)/den2

    return bij
            

def calcDen(N):

    den = 0.0
    for i in range(-N,N+1):  # i=-N,N
        den += (numpy.exp(-numpy.pi*abs( float(i)) / (float(N) / two ) ) )**2.0

    return numpy.sqrt(den)






    
    
    
