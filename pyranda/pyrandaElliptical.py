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
import scipy.sparse
from scipy.sparse.linalg import factorized,bicgstab,cg


class pyrandaPoisson:

    def __init__(self,BCtype,pysim,method="direct",itermax=100):
        """
        Solve an elliptical Poisson equation
        \Delta \phi = f(x)
        Second order FD on structured mesh assumed
        """
        self.BCtype = BCtype
        self.pysim = pysim
        self.method = method
        self.iterMax = itermax
        
        nx = self.nx = pysim.mesh.options['nn'][0]
        ny = self.ny = pysim.mesh.options['nn'][1]
        nz = self.nz = pysim.mesh.options['nn'][2]
        
        dx = (pysim.mesh.options['xn'][0] - pysim.mesh.options['x1'][0])/nx

        self.dx = dx
        
        n = self.n = nx*ny
        d = numpy.ones(n)
        b = numpy.zeros(n)

        d0 = d.copy()*-4.
        d1_lower = d.copy()[0:-1]
        d1_upper = d1_lower.copy()

        dnx_lower = d.copy()[0:-nx]
        dnx_upper = dnx_lower.copy()
            
        # Try to set both sides as Nuemann
        d1_upper[nx-1::nx] = 0.
        d1_upper[::nx]     = 2.
        
        d1_lower[nx-1::nx] = 0.
        d1_lower[nx-2::nx] = 2.
    
        dnx_upper[0:nx]    = 2.
        dnx_lower[-nx:]    = 2.

        # Avoid singular point
        d0[nx+1] = 1.0
        d1_upper[nx+1] = 0.0
        dnx_upper[nx+1] = 0.0
        d1_lower[nx] = 0.0
        dnx_lower[1] = 0.0
        
        d0 /= (dx*dx)
        d1_upper /= (dx*dx)
        d1_lower /= (dx*dx)
        dnx_upper /= (dx*dx)
        dnx_lower /= (dx*dx)

        self.solver=None
        self.A = None
        
        if method == 'direct':

            if (self.nz > 1 ):
                pysim.iprint("Error: Only 2D problems supported")
                exit()
            
            A = scipy.sparse.diags([d0, d1_upper, d1_lower, dnx_upper, dnx_lower], [0, 1, -1, nx, -nx], format='csc')
            
            self.solver = factorized(A)
            self.A = A
            

    def solve(self,rhs,guess):


        if self.method == 'direct':
            
            b = numpy.zeros(self.n) #RHS
            for j in range(0,self.ny):
                for i in range(0, self.nx):
                    b[j + i*self.ny] = rhs[i,j,0]


            sol = self.solver( b )
            mysol = rhs * 0.0

            for j in range(0,self.ny):
                for i in range(0, self.nx):
                    mysol[i,j,0] = sol[j + i*self.ny]

            return mysol

            
        if self.method == 'gauss':


            res = 1.0
            resMax = .01

            dx2 = self.dx * self.dx
            dy2 = 1.
            dz2 = 1.
            sixth = 1.0/6.0
            fourth = 1.0/4.0

            inter = numpy.index_exp[1:-1,1:-1,:]
            iF = numpy.index_exp[:-2,1:-1,:]
            iL = numpy.index_exp[2:,1:-1,:]
            
            jF = numpy.index_exp[1:-1,:-2,:]
            jL = numpy.index_exp[1:-1,2:,:]
            
            #inter = numpy.index_exp[1:-1,1:-1,1:-1]
            #iF = numpy.index_exp[:-2,1:-1,1:-1]
            #iL = numpy.index_exp[2:,1:-1,1:-1]
            
            #jF = numpy.index_exp[1:-1,:-2,1:-1]
            #jL = numpy.index_exp[1:-1,2:,1:-1]

            #kF = numpy.index_exp[1:-1,1:-1,:-2]
            #kL = numpy.index_exp[1:-1,1:-1,2:]

            #mysol[inter] = guess
            mysol = self.pysim.PyMPI.ghost( guess , np=1 , clip=False)
            
            old = mysol*1.0
            cnt = 1

            periodicX = self.pysim.PyMPI.periodic[0]
            periodicY = self.pysim.PyMPI.periodic[1]
            periodicZ = self.pysim.PyMPI.periodic[2]
            
            # Do stupid gauss siedel relaxation method
            while abs(res) > resMax:

                mysol[inter] = ( (mysol[iF] + mysol[iL]) +
                                 (mysol[jF] + mysol[jL])
                                 #(mysol[kF] + mysol[kL]) -
                                 - rhs * dx2 ) * fourth


                
                # Comm ghost data
                mysol = self.pysim.PyMPI.ghost( mysol[inter], np=1 , clip=False)

                if not periodicX:
                    if self.pysim.PyMPI.x1proc:
                        mysol[0,:,:] = mysol[2,:,:]
                    if self.pysim.PyMPI.xnproc:
                        mysol[-1,:,:] = mysol[-3,:,:]

                if not periodicY:
                    if self.pysim.PyMPI.y1proc:
                        mysol[:,0,:] = mysol[:,2,:]
                    if self.pysim.PyMPI.ynproc:
                        mysol[:,-1,:] = mysol[:,-3,:]
                
                cnt += 1
                if ( cnt > self.iterMax):
                    break

            return mysol[inter]
                
                
