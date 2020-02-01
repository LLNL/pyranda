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

    def __init__(self,BCtype,pysim):
        """
        Solve an elliptical Poisson equation
        \Delta \phi = f(x)
        Second order FD on structured mesh assumed
        """
        self.BCtype = BCtype
        self.pysim = pysim
        
        nx = self.nx = pysim.mesh.options['nn'][0]
        ny = self.ny = pysim.mesh.options['nn'][1]
        nz = self.nz = pysim.mesh.options['nn'][2]

        if (self.nz > 1 ):
            pysim.iprint("Error: Only 2D problems supported")
        
        dx = (pysim.mesh.options['xn'][0] - pysim.mesh.options['x1'][0])/nx 
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
        #d0[nx+1] = 1.0
        #d1_upper[n/2] = 0.0
        #dnx_upper[n/2] = 0.0
        #d1_lower[nx] = 0.0
        #dnx_lower[1] = 0.0
        
        d0 /= (dx*dx)
        d1_upper /= (dx*dx)
        d1_lower /= (dx*dx)
        dnx_upper /= (dx*dx)
        dnx_lower /= (dx*dx)
        A = scipy.sparse.diags([d0, d1_upper, d1_lower, dnx_upper, dnx_lower], [0, 1, -1, nx, -nx], format='csc')

        self.solver = factorized(A)
        self.A = A
        

    def solve(self,rhs):

        b = numpy.zeros(self.n) #RHS
        for j in range(0,self.ny):
            for i in range(0, self.nx):
                b[j + i*self.ny] = rhs[i,j,0]

        #b[0] = 0.0
        sol = self.solver( b )
        #sol = bicgstab(self.A, b )[0]
        #sol = cg(self.A, b )[0]
        mysol = rhs * 0.0

        #import pdb
        #pdb.set_trace()
        
        for j in range(0,self.ny):
            for i in range(0, self.nx):
                mysol[i,j,0] = sol[j + i*self.ny]

        mysol -= mysol.mean()
        return mysol
                
                
