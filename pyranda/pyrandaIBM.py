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

immersed_iter = 2
immersed_CFL = 0.5
immersed_EPS = 0.1

class pyrandaIBM(pyrandaPackage):

    def __init__(self,pysim):

        PackageName = 'IBM'
        pyrandaPackage.__init__(self,PackageName,pysim)


    def get_sMap(self):
        sMap = {}
        sMap['ibmWall('] = "self.packages['IBM'].ibmWall("
        sMap['ibmV('] = "self.packages['IBM'].ibmVel("
        sMap['ibmS('] = "self.packages['IBM'].ibmS("
        self.sMap = sMap

        
                 
   def ibmWallModel(self,vel,phi,gphi,global_u0,phivar=None):
       # do we have dx == dy == dz ????
       dx = self.pyranda.mesh.dx
       lo = self.pyranda.PyMPI.chunk_3d_lo
       ax, ay, az = phi.shape
       u0 = numpy.zeros((ax, ay, az))

       phix = gphi[0]
       phiy = gphi[1]
       phiz = gphi[2]

       tmp0 = (h_wm-phi)*phix
       tmp1 = (h_wm-phi)*phiy
       tmp2 = (h_wm-phi)*phiz

       shift0 = (tmp0/dx).astype(int)
       shift1 = (tmp1/dy).astype(int)
       shift2 = (tmp2/dz).astype(int)

       offset0 = (tmp0 - dx*shift0)/dx
       offset1 = (tmp1 - dy*shift1)/dy
       offset2 = (tmp2 - dz*shift2)/dz

       def interp(C):
           c00 = C[0,0,0]*(1-offset0)+C[1,0,0]*offset0
           c01 = C[0,0,1]*(1-offset0)+C[1,0,1]*offset0
           c10 = C[0,1,0]*(1-offset0)+C[1,1,0]*offset0
           c11 = C[0,1,1]*(1-offset0)+C[1,1,1]*offset0

           c0 = c00*(1-offset1)+c10*offset1
           c1 = c01*(1-offset1)+c11*offset1

           c = c0*(1-offset2)+c1*offset2
           return c

       for i in range(ax):
           for j in range(ay):
               for k in range(az):
                   ii = lo[0]+i+shift0
                   jj = lo[1]+j+shift1
                   kk = lo[2]+k+shift2
                   neighbs = global_u0[ii:ii+2,jj:jj+2,kk:kk+2]
                   u0[i,j,k] = interp(neighbs)

       # implement slip velocity and enforce wall model bc
       # essentially copying the slip velocity function with slip=True
       v1 = vel[0]*1.0
       v2 = vel[1]*1.0
       v3 = vel[2]*1.0

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
       v3 = self.smooth_terrain(SDF,gDx,gDy,gDz,v3,0.0,new=new)

       normal = v1*gDx+v2*gDy+v3*gDz

       vn =  numpy.where( SDF < lens, normal, 0.0 )

       # Remove normal velocity
       v1 = v1 - vn*gDx
       v2 = v2 - vn*gDy
       v3 = v3 - vn*gDz

       # Normalize the tangent vector
       norm = numpy.sqrt( v1*v1+v2*v2+v3*v3 )
       u1 = v1/norm
       u2 = v2/norm
       u3 = v3/norm

       # set velocity to zero for phi < 0
       interior_mask = numpy.where( SDF < lens, 0.0, 1.0)
       v1 = v1*interior_mask
       v2 = v2*interior_mask
       v3 = v3*interior_mask

       # Compute linear velocity through zero level
       tmp = numpy.where( SDF < lens, 0.0 , normal/(SDF) )
       tmp_n = self.smooth_terrain(SDF,gDx,gDy,gDz,tmp,lens,new=False)
       vn = numpy.where( SDF < lens, tmp_n*SDF, )

       # Compute the linear velocity through the wall model bc
       tmp = numpy.where( SDF < lens, u0 , normal/(SDF) )
       tmp_t = self.smooth_terrain(SDF,gDx,gDy,gDz,tmp,lens,new=False)
       vt = numpy.where( SDF < lens, tmp_t*SDF, norm )

       # Add velocity linear profile in tangential and normal directions
       v1 = v1 + vt*u1 + vn*gDx
       v2 = v2 + vt*u2 + vn*gDy
       v3 = v3 + vt*u3 + vn*gDz

       if phivar:
           v1 += v1_phi
           v2 += v2_phi
           v3 += v3_phi

       return [v1, v2, v3]

    def ibmVel(self,vel,phi,gphi,phivar=None):

        u = vel[0]
        v = vel[1]
        w = vel[2]
    
        phix = gphi[0]
        phiy = gphi[1]
        phiz = gphi[2]

        lens = self.pyranda.mesh.GridLen * immersed_EPS
    
        return self.slip_velocity( phi ,phix,phiy,phiz,
                                   u,v,w,lens,new=False,phivar=phivar,slip=True)

    def ibmWall(self,vel,phi,gphi,phivar=None):

        u = vel[0]
        v = vel[1]
        w = vel[2]
    
        phix = gphi[0]
        phiy = gphi[1]
        phiz = gphi[2]

        lens = self.pyranda.mesh.GridLen * immersed_EPS
    
        return self.slip_velocity( phi ,phix,phiy,phiz,
                                   u,v,w,lens,phivar=None,slip=False)
        
    
    def ibmS(self,scalar,phi,gphi):
        
        phix = gphi[0]
        phiy = gphi[1]
        phiz = gphi[2]
        
        epsi = 0.0    
    
        return self.smooth_terrain( phi, phix, phiy, phiz,
                                    scalar,epsi)
    
        

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


    def slip_velocity(self,SDF,gDx,gDy,gDz,v1_in,v2_in,v3_in,lens,new=False,phivar=None,slip=True):

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
        v3 = self.smooth_terrain(SDF,gDx,gDy,gDz,v3,0.0,new=new)
            

                
        if slip:


            norm = v1*gDx+v2*gDy+v3*gDz
            
            vn =  numpy.where( SDF < lens, norm, 0.0 )
            tmp = numpy.where( SDF < lens, 0.0 , norm/(SDF) )            
            
            # Remove normal velocity
            v1 = v1 - vn*gDx
            v2 = v2 - vn*gDy
            v3 = v3 - vn*gDz
            
            # Compute linear velocity through zero level
            tmp = self.smooth_terrain(SDF,gDx,gDy,gDz,tmp,lens,new=new)
            vn = numpy.where( SDF < lens, tmp*SDF, 0.0 )
        
            # Add velocity linear profile
            v1 = v1 + vn*gDx
            v2 = v2 + vn*gDy
            v3 = v3 + vn*gDz

            # Based on reconstructing a velocity gradient across the interface
            #dudn = Cf Re U / L  #  ~ \tau_w / mu  ... Cf = .004, L = Re scale, u\infty
            #dudn = .008 * 1.0e5 * 4.0 / 2.0             
            #norm   = numpy.sqrt( v1*v1+v2*v2+v3*v3 ) + 1.0e-10
            #factor = (1.0 + SDF*dudn/norm)  # 0.9
            if False:

                # Reduce magnitude inside... keep vector
                # ... This is wrong and should scale with SDF somehow...
                v1 = numpy.where( SDF < lens, v1*.9, v1)
                v2 = numpy.where( SDF < lens, v2*.9, v2)
                v3 = numpy.where( SDF < lens, v3*.9, v3)
                
                factor  = (1.0 + SDF*.0001)
            
                v1 = numpy.where( SDF < lens, v1*factor, v1)
                v2 = numpy.where( SDF < lens, v2*factor, v2)
                v3 = numpy.where( SDF < lens, v3*factor, v3)


            
            
            

            
        else:
            norm = numpy.sqrt( v1*v1+v2*v2+v3*v3 )
            
            vn =  numpy.where( SDF < lens, norm, 0.0 )
            tmp = numpy.where( SDF < lens, 0.0 , norm/(SDF) )            

            # temp vectors
            inorm = 1.0 / (norm + 1.0e-16)
            tmpx = v1 * inorm
            tmpy = v2 * inorm
            tmpz = v3 * inorm
            
            v1 -= vn*tmpx
            v2 -= vn*tmpy
            v3 -= vn*tmpz

            # Compute linear velocity through zero level
            tmp = self.smooth_terrain(SDF,gDx,gDy,gDz,tmp,lens,new=new)
            vn = numpy.where( SDF < lens, tmp*SDF, 0.0 )
            
            # Add velocity linear profile
            v1 += vn*tmpx
            v2 += vn*tmpy
            v3 += vn*tmpz
            
            
        if phivar:
            v1 += v1_phi
            v2 += v2_phi
            v3 += v3_phi

        return [v1,v2,v3]

    
