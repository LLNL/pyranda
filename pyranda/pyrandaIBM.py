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

immersed_iter = 5
immersed_CFL = 0.5
immersed_EPS = 0.5

class pyrandaIBM(pyrandaPackage):

    def __init__(self,pysim):

        PackageName = 'IBM'
        pyrandaPackage.__init__(self,PackageName,pysim)


    def get_sMap(self):
        sMap = {}
        sMap['ibmV('] = "self.packages['IBM'].ibmVel("
        sMap['ibmS('] = "self.packages['IBM'].ibmS("
        sMap['ibmS2('] = "self.packages['IBM'].ibmS2("
        sMap['ibmC('] = "self.packages['IBM'].ibmC("
        self.sMap = sMap
        
                 
    def ibmVel(self,vel,phi,gphi,phivar=None):

        u = vel[0]
        v = vel[1]
        w = vel[2]
    
        phix = gphi[0]
        phiy = gphi[1]
        phiz = gphi[2]

        lens = self.pyranda.mesh.GridLen * immersed_EPS
    
        return self.slip_velocity( phi ,phix,phiy,phiz,
                                   u,v,w,lens,new=False,phivar=phivar)
    
    def ibmS(self,scalar,phi,gphi):
        
        phix = gphi[0]
        phiy = gphi[1]
        phiz = gphi[2]
        
        epsi = 0.0    
    
        return self.smooth_terrain( phi, phix, phiy, phiz,
                                    scalar,epsi)


    def ibmC(self,scalar,SDF,gphi,Cval):

        SS = scalar - Cval
        lens = self.pyranda.mesh.GridLen * immersed_EPS
        gDx = gphi[0]
        gDy = gphi[1]
        gDz = gphi[2]
        
        #vn =  numpy.where( SDF < lens, SS, 0.0 )
        tmp = numpy.where( SDF < lens, 0.0 , SS/SDF )            
        tmp = numpy.where( SDF > -lens, tmp , SS/SDF )
        #tmp = numpy.where( SDF == 0.0, 0.0,SS/SDF)
        
        
        # Compute linear velocity through zero level
        #lens = 0.0
        tmp = self.smooth_terrain(SDF,gDx,gDy,gDz,tmp,lens)
        vn = numpy.where( SDF < lens, tmp*SDF + Cval, scalar )

        return vn
    
    
    def ibmS2(self,scalar,phi,gphi):
        
        phix = gphi[0]
        phiy = gphi[1]
        phiz = gphi[2]
        
        epsi = 0.0    

        # Commpute dV/dPhi
        [tvx,tvy,tvz] = self.pyranda.grad(scalar)
        dvdp =  tvx * phix #+ tvy*phiy + tvz*phiz
        dvdp = numpy.where( phi > 0.0 , dvdp, 0.0 )

        # This gradient goes to zero
        Sdvdp = self.smooth_terrain( phi, phix,phiy,phiz,dvdp,epsi)
        #Sdvdp = self.ibmC( dvdp, phi, gphi, 0.0)

        # Only apply der near zero level
        #GridLen = self.pyranda.mesh.GridLen
        #Sdvdp *= numpy.where( numpy.abs( phi ) < 5.0*GridLen, 1.0, 0.0 )
        
        Sv = self.smooth_terrain( phi, phix, phiy, phiz,
                                  scalar,epsi)

        
        S = numpy.where( phi <= epsi, Sv + phi * Sdvdp, Sv)
        return S
        
        

    def smooth_terrain(self,SDF,gDx,gDy,gDz,val_in,epsi,new=False):
        
        val = val_in * 1.0
        GridLen = self.pyranda.mesh.GridLen
        
        for i in range(immersed_iter):
            [tvx,tvy,tvz] = self.pyranda.grad(val)
            term = tvx*gDx+tvy*gDy+tvz*gDz
            #term += self.pyranda.laplacian(SDF)*val
            val = numpy.where( SDF <= epsi , val + immersed_CFL*GridLen*term , val )
            Tval = self.pyranda.gfilter(val)
            val = numpy.where( SDF <= epsi , Tval, val )
        
        return val


    def smooth_terrain2(self,SDF,gDx,gDy,gDz,val_in,epsi,new=False):
        
        val = val_in * 1.0
        GridLen = self.pyranda.mesh.GridLen
        
        for i in range(immersed_iter):
            [tvx,tvy,tvz] = self.pyranda.grad(val)
            term = tvx*gDx+tvy*gDy+tvz*gDz
            val = numpy.where( SDF <= epsi , val + immersed_CFL*GridLen*term , val )
            Tval = self.pyranda.gfilter(val)
            val = numpy.where( SDF <= epsi , Tval, val )
        
        return val

    def slip_velocity(self,SDF,gDx,gDy,gDz,v1_in,v2_in,v3_in,lens,new=False,phivar=None):

        v1 = v1_in*1.0
        v2 = v2_in*1.0
        v3 = v3_in*1.0

        if phivar:
            v1_phi = phivar[0]
            v2_phi = phivar[1]
            v3_phi = phivar[2]

            # Transform to interface velocity
            v1 -= v1_phi
            v2 -= v2_phi
            v3 -= v3_phi

        v1 = self.smooth_terrain(SDF,gDx,gDy,gDz,v1,0.0,new=new)
        v2 = self.smooth_terrain(SDF,gDx,gDy,gDz,v2,0.0,new=new)
        #v3 = self.smooth_terrain(SDF,gDx,gDy,gDz,v3,0.0,new=new)
            
        norm = v1*gDx+v2*gDy+v3*gDy
        vn =  norm #numpy.where( SDF < lens, norm, 0.0 )
    
        # Remove normal velocity
        v1 = v1 - vn*gDx
        v2 = v2 - vn*gDy
        v3 = v3 - vn*gDz

        # Compute linear velocity through zero level
        #tmp = numpy.where( SDF < lens, 0.0 , norm/SDF )
        #tmp = numpy.where( SDF > -lens, tmp , norm/SDF )
        #tmp = self.smooth_terrain(SDF,gDx,gDy,gDz,tmp,lens,new=new)
        #vn = numpy.where( SDF < lens, tmp*SDF, 0.0 )
        vn = self.ibmC(norm,SDF,[gDx,gDy,gDz],0.0)
    
        # Add velocity linear profile
        v1 = v1 + vn*gDx
        v2 = v2 + vn*gDy
        v3 = v3 + vn*gDz

        if phivar:
            v1 += v1_phi
            v2 += v2_phi
            v3 += v3_phi

        return [v1,v2,0.0]

    
