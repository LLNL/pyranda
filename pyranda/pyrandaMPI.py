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
from functools import reduce

from . import parcop

# Global patch/levels to create new instances
GLOBAL_PATCH = 0
GLOBAL_LEVEL = 0
LAST_PATCH = -1
LAST_LEVEL = 0
PATCH_MAX = 100
LEVEL_MAX = 10

class pyrandaMPI():

    def __init__(self,mesh):


        meshOptions = mesh.options
        self.nx = meshOptions['nn'][0]
        self.ny = meshOptions['nn'][1]
        self.nz = meshOptions['nn'][2]
        
        x1 = meshOptions['x1'][0]
        xn = meshOptions['xn'][0]
        y1 = meshOptions['x1'][1]
        yn = meshOptions['xn'][1]
        z1 = meshOptions['x1'][2]
        zn = meshOptions['xn'][2]
        
        dx = (meshOptions['xn'][0]-meshOptions['x1'][0])/max(self.nx-1,1)
        dy = (meshOptions['xn'][1]-meshOptions['x1'][1])/max(self.ny-1,1)
        dz = (meshOptions['xn'][2]-meshOptions['x1'][2])/max(self.nz-1,1)
        self.dx = dx
        self.dy = dy
        self.dz = dz
        
        periodic = meshOptions['periodic']
        self.periodic = periodic

        self.coordsys = meshOptions['coordsys']

        # Get the mesh function if it exists
        try:
            self.meshFunction = meshOptions['meshFunction']
        except:
            self.meshFunction = None
            
        self.comm  = MPI.COMM_WORLD
        self.fcomm = self.comm.py2f()
            
            
        self.order = (10,10,10)
        self.filter_type = ('compact', 'compact', 'compact')
        

        # Compute the domain decomposition
        [px,py,pz] = decomp( self.comm.Get_size(),
                             self.nx, self.ny, self.nz )

        if (px*py*pz != self.comm.Get_size()):
            if self.comm.Get_rank() == 0:
                print("Tried (px,py,pz) = (%s,%s,%s)" % (px,py,pz))
                print("Nprocs = %s" % (px*py*pz))
                print("Available procs = %s" % self.comm.Get_size())
                print('Error: Processor under/over utilization not allowed')
                exit()
                
        
        
        bx1 = "NONE"
        bxn = "NONE"
        by1 = "NONE"#"PERI"
        byn = "NONE"#"PERI"
        bz1 = "NONE"#"PERI"
        bzn = "NONE"#"PERI"
        if periodic[0]:
            bx1 = "PERI"
            bxn = "PERI"
        if periodic[1]:
            by1 = "PERI"
            byn = "PERI"
        if periodic[2]:
            bz1 = "PERI"
            bzn = "PERI"
        
        # Set at unique patch/level for each instance
        global GLOBAL_PATCH
        global LAST_PATCH
        global GLOBAL_LEVEL
        global LAST_LEVEL
        
        GLOBAL_PATCH = LAST_PATCH + 1
        if GLOBAL_PATCH >= PATCH_MAX:
            LAST_LEVEL = LAST_LEVEL + 1
            GLOBAL_PATCH = 0
            
        GLOBAL_LEVEL = int(LAST_LEVEL)
        self.patch = int(GLOBAL_PATCH)
        self.level = int(GLOBAL_LEVEL)

        parcop.parcop.setup( self.patch, self.level , self.fcomm,
                             self.nx,self.ny,self.nz, 
                             px,py,pz,self.coordsys,
                             x1,xn,y1,yn,z1,zn,
                             bx1,bxn,by1,byn,bz1,bzn)


        
        self.chunk_3d_size = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_lo   = numpy.zeros(3, dtype=numpy.int32, order='F')
        self.chunk_3d_hi   = numpy.zeros(3, dtype=numpy.int32, order='F')

        self.setPatch()
        self.xcom  = MPI.Comm.f2py( parcop.parcop.commx()  )
        self.ycom  = MPI.Comm.f2py( parcop.parcop.commy()  )
        self.zcom  = MPI.Comm.f2py( parcop.parcop.commz()  )
        self.xycom = MPI.Comm.f2py( parcop.parcop.commxy() )
        self.yzcom = MPI.Comm.f2py( parcop.parcop.commyz() )
        self.xzcom = MPI.Comm.f2py( parcop.parcop.commxz() )

            

        self.chunk_3d_size[0] = self.ax = int( self.nx / px )
        self.chunk_3d_size[1] = self.ay = int( self.ny / py )
        self.chunk_3d_size[2] = self.az = int( self.nz / pz )

        self.chunk_3d_lo[0] = self.xcom.rank     * self.chunk_3d_size[0] + 1
        self.chunk_3d_hi[0] = (self.xcom.rank+1) * self.chunk_3d_size[0]

        self.chunk_3d_lo[1] = self.ycom.rank     * self.chunk_3d_size[1] + 1
        self.chunk_3d_hi[1] = (self.ycom.rank+1) * self.chunk_3d_size[1]
            
        self.chunk_3d_lo[2] = self.zcom.rank     * self.chunk_3d_size[2] + 1
        self.chunk_3d_hi[2] = (self.zcom.rank+1) * self.chunk_3d_size[2]

        self.chunk_3d_lo = self.chunk_3d_lo - 1 # Convert to 0 based indexing
        self.chunk_3d_hi = self.chunk_3d_hi - 1 # Convert to 0 based indexing


        self.der  = parcop_der()
        self.fil  = parcop_sfil()
        self.gfil = parcop_gfil()


        #self.setPatch()        

        

        # Communicators and proc boundaries
        self.master = False
        if self.comm.rank == 0:
            self.master = True

        self.x1proc = False
        if self.xcom.rank == 0:
            self.x1proc = True

        self.xnproc = False
        if self.xcom.rank == self.xcom.size - 1:
            self.xnproc = True

        self.y1proc = False
        if self.ycom.rank == 0:
            self.y1proc = True

        self.ynproc = False
        if self.ycom.rank == self.ycom.size - 1:
            self.ynproc = True
                                 
    def setPatch(self):                                 
        parcop.parcop.set_patch( self.patch, self.level )

    def emptyScalar(self):
        return numpy.zeros( self.chunk_3d_size, dtype=numpy.float64, order='F')*0.0

    def emptyVector(self):
        blk_size = numpy.append(self.chunk_3d_size,3)
        return numpy.zeros( blk_size, dtype=numpy.float64, order='F')*0.0

    
    def gather2D(self,idata,com,n1,n2,g1,g2):
        """
        Give idata(a1,a2) return the assembled
          global array of gdata(n1,n2)
        """        

        ldata = numpy.zeros( (n1,n2) )
        ldata[ g1[0]:g1[1], g2[0]:g2[1] ] = idata

        fldata = ldata.flatten()
        fldataT = fldata * 0.0
        com.Allreduce( fldata, fldataT, op=MPI.SUM )
        
        return numpy.reshape(fldataT, (n1,n2), order = 'C' )


    def sum2D(self,data,com,n2,n3,index,g2,g3):

        a2 = g2[1] - g2[0]
        a3 = g3[1] - g3[0]
        gsum = numpy.zeros((a2,a3))

        # Get the local proc mean
        lsum = numpy.sum( data, index )

        tsum = self.gather2D(lsum,com,n2,n3,g2,g3)

        return tsum

    def sum3D(self,data):

        lsum = numpy.sum( data )
        tsum = self.comm.allreduce( lsum, op=MPI.SUM )
        
        return tsum

    def max3D(self,data):

        lmax = numpy.max( data )
        tmax = self.comm.allreduce( lmax, op=MPI.MAX )

        return tmax

    def min3D(self,data):

        lmin = numpy.min( data )
        tmin = self.comm.allreduce( lmin, op=MPI.MIN )

        return tmin

    def xbar(self,data):        
        return self.sum2D( data, self.comm,self.ny,self.nz,0,
                           [ self.chunk_3d_lo[1] , self.chunk_3d_hi[1] ],
                           [ self.chunk_3d_lo[2] , self.chunk_3d_hi[2] ]  )
  
    def ybar(self,data):
        return self.sum2D( data, self.comm,self.nx,self.nz,1,
                           [self.chunk_3d_lo[0] , self.chunk_3d_hi[0]],
                           [self.chunk_3d_lo[2] , self.chunk_3d_hi[2]])
    
    def zbar(self,data):
        return self.sum2D( data, self.comm,self.nx,self.ny,2,
                           [self.chunk_3d_lo[0],self.chunk_3d_hi[0]+1 ],
                           [self.chunk_3d_lo[1],self.chunk_3d_hi[1]+1 ])
    




class parcop_der:

    def __init__(self):        
        pass

    def ddx(self,val):        
        return parcop.parcop.ddx( val )

    def ddy(self,val):        
        return parcop.parcop.ddy( val )

    def ddz(self,val):        
        return parcop.parcop.ddz(  val )

    def div(self,fx,fy,fz):
        return parcop.parcop.divergence(fx,fy,fz)

    def laplacian(self,val):        
        return parcop.parcop.plaplacian(  val )

    def ring(self,val):
        return parcop.parcop.pring(  val )

class parcop_mesh:

    def __init__(self):
        pass

    #def get
    
class parcop_gfil:

    def __init__(self):
        pass

    def filter(self,val):
        return parcop.parcop.gfilter( val)

class parcop_sfil:

    def __init__(self):
        pass

    def filter(self,val):
        return parcop.parcop.sfilter(  val)




def factors(n):    
    ff = set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(pow(n, 0.5) + 1)) if n % i == 0)))
    return list(ff)


class pyDecomp:
    def __init__(self,xp,yp,zp,np,nx,ny,nz):

        self.xp = xp
        self.yp = yp
        self.zp = zp
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.np = np

        self.totalProcs = self.xp * self.yp * self.zp

        ax = nx / xp
        ay = ny / yp
        az = nz / zp

        surface = 2*ax*ay + 2*ay*az + 2*az*ax
        volume  = ax*ay*az

        self.surf2vol = float(surface) / float(volume)

        # optimal surface to volume ratio
        sopt = float(volume)**(2./3.) * 6.0  / float( volume )
        

        self.score = 1.0
        
        # Too many procs
        if self.totalProcs > self.np:
            self.score *= 0.0            
            

        # Too small in any direction
        if (ax > 1) and (ax < 16):
            self.score *= 0.0

        if (ay > 1) and (ay < 16):
            self.score *= 0.0

        if (az > 1) and (az < 16):
            self.score *= 0.0
            

        # Total procs (under-utilization penalty
        #self.score *= self.totalProcs / self.np
        if self.totalProcs < self.np:
            self.score *= 0.0
        

        # Surface/Volume
        safrac = .1
        self.score -= safrac*(1.0 - sopt / self.surf2vol  )
        
            
        
        

def decomp(nprocs,nx,ny,nz):
    """ 
    Return the optimal domain decomp for nprocs
    """

    # Get factors in each direction
    xfac = factors(nx)
    yfac = factors(ny)
    zfac = factors(nz)
    

    # Get all combos
    decomps = []
    for xf in xfac:
        for yf in yfac:
            for zf in zfac:
                xproc = nx / xf
                yproc = ny / yf
                zproc = nz / zf
                decomps.append( pyDecomp( xproc,yproc,zproc,nprocs,nx,ny,nz) )
                
    decomps.sort(key=lambda x: x.score, reverse=True)    
    return [decomps[0].xp,decomps[0].yp,decomps[0].zp]
