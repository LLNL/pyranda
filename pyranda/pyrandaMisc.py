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
import scipy.interpolate as SI
import scipy     
            
        
class ipoint:
    """
    Interpolated object to rapidly get values
      Note: grids should be ghosted in the +i,j+ directions
    """
    def __init__(self,xpt,ypt,xgrid,ygrid):
        

        self.value = 0.0
        self.onProc = False
        
        self.x = xpt
        self.y = ypt

        # Bounding box for determining proc location
        onX = ( (self.x > xgrid.min() ) and  (self.x < xgrid.max() ) )
        onY = ( (self.y > ygrid.min() ) and  (self.y < ygrid.max() ) )        
        if onX and onY:
            self.onProc = True
        else:
            return
            
        self.ii1 = -1
        self.iin = -1

        self.jj1 = -1
        self.jjn = -1

        [self.nx,self.ny,self.nz] = xgrid.shape
               
        self.weight  = []

        x_y = numpy.dstack([xgrid.ravel(),ygrid.ravel()])[0]

        pts = numpy.zeros( (2,1) )
        pts[0] = self.x*1.0
        pts[1] = self.y*1.0
        points = list(pts.transpose())

        mytree = scipy.spatial.cKDTree(x_y)
        dist, indexes = mytree.query(points)

        ii = int(indexes[0]/self.ny)
        jj = indexes[0]%self.ny

        x0 = numpy.abs(xgrid[ii,jj] - self.x)
        y0 = numpy.abs(ygrid[ii,jj] - self.y)

        self.ii0 = ii
        self.jj0 = jj

        self.xlocal = numpy.array([ [xgrid[ii-1,jj-1],xgrid[ii,jj-1],xgrid[ii+1,jj-1]],
                                    [xgrid[ii-1,jj  ],xgrid[ii,jj  ],xgrid[ii+1,jj  ]],
                                    [xgrid[ii-1,jj+1],xgrid[ii,jj+1],xgrid[ii+1,jj+1]] ])

        self.ylocal = numpy.array([ [ygrid[ii-1,jj-1],ygrid[ii,jj-1],ygrid[ii+1,jj-1]],
                                    [ygrid[ii-1,jj  ],ygrid[ii,jj  ],ygrid[ii+1,jj  ]],
                                    [ygrid[ii-1,jj+1],ygrid[ii,jj+1],ygrid[ii+1,jj+1]] ])

        #self.xlocal = self.xlocal[1:,1:]
        #self.ylocal = self.ylocal[1:,1:]
        #self.xlocal = self.xlocal[:-1,:-1]
        #self.ylocal = self.ylocal[:-1,:-1]
        


        # Fast binlinear version
        indices = [ [ii,jj] ]
        indices.append( [ii-1,jj] )
        indices.append( [ii,jj-1] )
        indices.append( [ii-1,jj-1] )

        indices = [ [ii-1,jj-1],[ii,jj-1],[ii+1,jj-1],
                    [ii-1,jj  ],[ii,jj  ],[ii+1,jj  ],
                    [ii-1,jj+1],[ii,jj+1],[ii+1,jj+1] ]
        
        xD = [ ]
        yD = [ ]
        for ind in indices:
            xD.append( xgrid[ind[0],ind[1]] )
            yD.append( ygrid[ind[0],ind[1]] )


        D = numpy.sqrt( (numpy.array(xD)-self.x)**2 + (numpy.array(yD)-self.y)**2 )
        #D = D**2.0
        D = 1.0 / (D+1.0e-12)
        DT = numpy.sum( D )

        self.weights = []
        for Di in D:
            self.weights.append( Di / DT )

        self.indices = indices
            
            
        

    def interp(self,vals,imethod='linear'):


        if (not self.onProc):
            return 0.0
        
        ii = self.ii0
        jj = self.jj0
        vlocal = numpy.array([ [vals[ii-1,jj-1],vals[ii,jj-1],vals[ii+1,jj-1]],
                             [vals[ii-1,jj  ],vals[ii,jj  ],vals[ii+1,jj  ]],
                             [vals[ii-1,jj+1],vals[ii,jj+1],vals[ii+1,jj+1]] ])

        #vlocal = vlocal[1:,1:]
        #vlocal = vlocal[:-1,:-1]
        
        pts = [ numpy.array([self.x]), numpy.array([self.y]) ]
        dc = getPoints(self.xlocal,self.ylocal,vlocal,pts,imethod)

        return dc

    def interp1(self,vals,imethod='linear'):


        if (not self.onProc):
            return 0.0

        dc = 0.0
        for ind,wi in zip(self.indices,self.weights):
            dc += wi[0]*vals[ind[0],ind[1]]
            
        return dc

    
def getPoints(x,y,f,pts,imethod):

    xx = x.ravel()
    yy = y.ravel()
    ff = f.ravel()
    xpts = pts[0].ravel()
    ypts = pts[1].ravel()
    
    dc = SI.griddata(numpy.array([xx,yy]).T,ff,numpy.array([xpts,ypts]).T,method=imethod)
        
    return dc[0]
