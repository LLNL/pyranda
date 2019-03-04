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


class pyrandaBC(pyrandaPackage):

    def __init__(self,pysim):

        PackageName = 'BC'
        pyrandaPackage.__init__(self,PackageName,pysim)

        self.bcList = {}

        self.BCdata = {}

        

    def get_sMap(self):
        sMap = {}
        sMap['bc.extrap('] =  "self.packages['BC'].extrap("
        sMap['bc.const(']  =  "self.packages['BC'].const("
        sMap['bc.exit(']   =  "self.packages['BC'].exitbc("
        sMap['bc.slip(']   =  "self.packages['BC'].slipbc("
        sMap['bc.field(']  =  "self.packages['BC'].applyBC('field',"
        sMap['bc.farfield(']  =  "self.packages['BC'].farfieldbc("
        self.sMap = sMap
        

    def applyBC(self,func,var,direction,val=None):

        pFunc = eval('self.' + func)
        
        if type(var) != type([]):
            var = [var]
        
        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                if (type(val) != type(None)):
                    pFunc( v, d , val )
                else:
                    pFunc( v, d  )
        
        
    def extrap(self,var,direction,order=2):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.extrapolate( v , d ,order)
            
        

    def extrapolate(self,var,direction,order):
        # Direction switch
        bcvar = None
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                if order == 2:
                    self.pyranda.variables[var].data[0,:,:] = 2*self.pyranda.variables[var].data[1,:,:] - self.pyranda.variables[var].data[2,:,:] 
                else:
                    self.pyranda.variables[var].data[0,:,:] = self.pyranda.variables[var].data[1,:,:] 

        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                if order == 2:
                    self.pyranda.variables[var].data[-1,:,:] = 2*self.pyranda.variables[var].data[-2,:,:] - self.pyranda.variables[var].data[-3,:,:]
                else:
                    self.pyranda.variables[var].data[-1,:,:] = self.pyranda.variables[var].data[-2,:,:]

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                if order == 2:
                    self.pyranda.variables[var].data[:,0,:] = 2*self.pyranda.variables[var].data[:,1,:] - self.pyranda.variables[var].data[:,2,:]
                else:
                    self.pyranda.variables[var].data[:,0,:] = self.pyranda.variables[var].data[:,1,:]

        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                if order == 2:
                    self.pyranda.variables[var].data[:,-1,:] = 2*self.pyranda.variables[var].data[:,-2,:] - self.pyranda.variables[var].data[:,-3,:]
                else:
                    self.pyranda.variables[var].data[:,-1,:] = self.pyranda.variables[var].data[:,-2,:]

                
    def const(self,var,direction,val):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.constant( v , d , val)
                
                
    def constant(self,var,direction,val):
        
        # Direction switch
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                self.pyranda.variables[var].data[0,:,:] = val
                
            
        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                self.pyranda.variables[var].data[-1,:,:] = val
            

        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                self.pyranda.variables[var].data[:,0,:] = val
                
            
        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                self.pyranda.variables[var].data[:,-1,:] = val

    def field(self,var,direction,field):

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

        

                 

    def slipbc(self,var,direction):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.slipvelbc( v , d )
                
    def slipvelbc(self,velocity,direction):
        """
        Allow tangential velocities
        """

        ax = self.pyranda.PyMPI.ax
        ay = self.pyranda.PyMPI.ay
        az = self.pyranda.PyMPI.az
        
        if direction == 'y1':

            if not self.BCdata.has_key('slipbc-y1'):
                
                x = self.pyranda.mesh.coords[0].data[:,:,:]
                y = self.pyranda.mesh.coords[1].data[:,:,:]
                z = self.pyranda.mesh.coords[2].data[:,:,:]

                #CASE(1) ! A direction
                dAdx = self.pyranda.getVar("dAx")
                dAdy = self.pyranda.getVar("dAy")
                xTan = dAdx/numpy.sqrt(dAdx**2+dAdy**2 )
                yTan = dAdy/numpy.sqrt(dAdx**2+dAdy**2 )

                # Trivial rotate 90 in 2D
                xN = -yTan
                yN =  xTan

                #CASE(2) ! B direction
                #xnormal = dBdx/numpy.sqrt(dBdx**2+dBdy**2 )#*sign
                #ynormal = dBdy/numpy.sqrt(dBdx**2+dBdy**2 )#*sign
                #znormal = dBdz/numpy.sqrt(dBdx**2+dBdy**2 )#*sign
                #CASE(3) ! C direction
                #xnormal = dCdx/numpy.sqrt(dCdx**2+dCdy**2+dCdz**2)*sign
                #ynormal = dCdy/numpy.sqrt(dCdx**2+dCdy**2+dCdz**2)*sign
                #znormal = dCdz/numpy.sqrt(dCdx**2+dCdy**2+dCdz**2)*sign
                
                norms = 0.0
                if self.pyranda.PyMPI.y1proc:

                    xb = x[:,1,:]
                    yb = y[:,1,:]
                    zb = z[:,1,:]

                    # Norms are assummed to be built into mesh!
                    norms = []
                    #mag = 0.0 * xb
                    #if ax > 1:
                    #    xbp = x[:,2,:]
                    #    norms.append( xbp - xb )
                    #    mag += (xbp - xb)**2
                    #if ay > 1:
                    #    ybp = y[:,2,:]
                    #    norms.append( ybp - yb )
                    #    mag += (ybp - yb)**2
                    #if az > 1:
                    #    zbp = z[:,2,:]
                    #    norms.append( zbp - zb )
                    #    mag += (zbp - zb)**2

                    #mag = numpy.sqrt( mag )
                    #for nn in range(len(norms)):
                    #    norms[nn] = norms[nn] / mag

                    norms.append( xN[:,0,:] )
                    norms.append( yN[:,0,:] )
                    
                        

                self.BCdata['slipbc-y1'] = norms
            else:
                norms = self.BCdata['slipbc-y1']

                
            # For free slip, always extrapolate first
            for uu in velocity:
                self.extrapolate(uu,'y1',order=2)
                
            if self.pyranda.PyMPI.y1proc:

                udotn = 0.0
                for uu in range(len(velocity)):
                    u = self.pyranda.variables[velocity[uu]].data[:,0,:]
                    udotn = udotn + u*norms[uu]

                # magnitude before
                mag0 = 0.0
                for uu in range(len(velocity)):
                    mag0 += self.pyranda.variables[velocity[uu]].data[:,0,:]**2
                mag0 = numpy.sqrt( mag0 ) #+ 1.0e-12
                    
                for uu in range(len(velocity)):                    
                    self.pyranda.variables[velocity[uu]].data[:,0,:] -= udotn*norms[uu]

                # magnitude after
                magF = 0.0
                for uu in range(len(velocity)):
                    magF += self.pyranda.variables[velocity[uu]].data[:,0,:]**2
                magF = numpy.sqrt( magF ) #+ 1.0e-12

                # Make sure magnitude is NOT larger
                for uu in range(len(velocity)):                    
                    self.pyranda.variables[velocity[uu]].data[:,0,:] = numpy.where(
                        magF > mag0,
                        self.pyranda.variables[velocity[uu]].data[:,0,:]*mag0/magF,
                        self.pyranda.variables[velocity[uu]].data[:,0,:] )

                
                

                    
    def exitbc(self,var,direction,norm=False):

        if type(var) != type([]):
            var = [var]

        if type(direction) != type([]):
            direction = [direction]

        for d in direction:
            for v in var:
                self.exit_boundary( v , d , norm)

    def exit_boundary(self,var,direction,norm):


        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                u1 = self.pyranda.variables[var].data[0,:,:]
                u2 = self.pyranda.variables[var].data[1,:,:]
                u3 = self.pyranda.variables[var].data[2,:,:]                
                self.pyranda.variables[var].data[0,:,:] = self.BENO(u1,u2,u3,norm)
        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                u1 = self.pyranda.variables[var].data[-1,:,:]
                u2 = self.pyranda.variables[var].data[-2,:,:]
                u3 = self.pyranda.variables[var].data[-3,:,:]                
                self.pyranda.variables[var].data[-1,:,:] = self.BENO(u1,u2,u3,norm)



    def BENO(self,u1,u2,u3,norm):
        d1 =        u2     - u1
        d2 = (2.0*u2-u3) - u1

        beno = numpy.where( numpy.abs(d1) <= numpy.abs(d2), u1+d1, u1+d2)
        beno = numpy.where(d1*d2 <= 0.0, u1, beno)
              
        if norm:        
            beno = numpy.maximum( numpy.minimum( beno, numpy.maximum(0.0,u2) ), numpy.minimum( 0.0, u2) )

        return beno



    def getnormal(self,direction):
        
        dAdx = self.pyranda.getVar("dAx")
        dAdy = self.pyranda.getVar("dAy")
        dAdz = self.pyranda.getVar("dAz")

        dBdx = self.pyranda.getVar("dBx")
        dBdy = self.pyranda.getVar("dBy")
        dBdz = self.pyranda.getVar("dBz")

        dCdx = self.pyranda.getVar("dCx")
        dCdy = self.pyranda.getVar("dCy")
        dCdz = self.pyranda.getVar("dCz")
        
        dataKey = "farfield-%s" % direction
        if self.BCdata.has_key( dataKey ):
            [xnormal,ynormal,znormal] = self.BCdata[dataKey]

        else:
            if (direction == 'x1'):
                xnormal = -dAdx[0,:,:]/numpy.sqrt(dAdx[0,:,:]**2+dAdy[0,:,:]**2+dAdz[0,:,:]**2)
                ynormal = -dAdy[0,:,:]/numpy.sqrt(dBdx[0,:,:]**2+dBdy[0,:,:]**2+dBdz[0,:,:]**2)
                znormal = -dAdz[0,:,:]/numpy.sqrt(dCdx[0,:,:]**2+dCdy[0,:,:]**2+dCdz[0,:,:]**2)
            if (direction == 'y1'):
                xnormal = -dBdx[:,0,:]/numpy.sqrt(dAdx[:,0,:]**2+dAdy[:,0,:]**2+dAdz[:,0,:]**2)
                ynormal = -dBdy[:,0,:]/numpy.sqrt(dBdx[:,0,:]**2+dBdy[:,0,:]**2+dBdz[:,0,:]**2)
                znormal = -dBdz[:,0,:]/numpy.sqrt(dCdx[:,0,:]**2+dCdy[:,0,:]**2+dCdz[:,0,:]**2)
            if (direction == 'z1'):
                xnormal = -dCdx[:,:,0]/numpy.sqrt(dAdx[:,:,0]**2+dAdy[:,:,0]**2+dAdz[:,:,0]**2)
                ynormal = -dCdy[:,:,0]/numpy.sqrt(dBdx[:,:,0]**2+dBdy[:,:,0]**2+dBdz[:,:,0]**2)
                znormal = -dCdz[:,:,0]/numpy.sqrt(dCdx[:,:,0]**2+dCdy[:,:,0]**2+dCdz[:,:,0]**2)

            if (direction == 'xn'):
                xnormal = dAdx[-1,:,:]/numpy.sqrt(dAdx[-1,:,:]**2+dAdy[-1,:,:]**2+dAdz[-1,:,:]**2)
                ynormal = dAdy[-1,:,:]/numpy.sqrt(dBdx[-1,:,:]**2+dBdy[-1,:,:]**2+dBdz[-1,:,:]**2)
                znormal = dAdz[-1,:,:]/numpy.sqrt(dCdx[-1,:,:]**2+dCdy[-1,:,:]**2+dCdz[-1,:,:]**2)
            if (direction == 'yn'):
                xnormal = dBdx[:,-1,:]/numpy.sqrt(dAdx[:,-1,:]**2+dAdy[:,-1,:]**2+dAdz[:,-1,:]**2)
                ynormal = dBdy[:,-1,:]/numpy.sqrt(dBdx[:,-1,:]**2+dBdy[:,-1,:]**2+dBdz[:,-1,:]**2)
                znormal = dBdz[:,-1,:]/numpy.sqrt(dCdx[:,-1,:]**2+dCdy[:,-1,:]**2+dCdz[:,-1,:]**2)
            if (direction == 'zn'):
                xnormal = dCdx[:,:,-1]/numpy.sqrt(dAdx[:,:,-1]**2+dAdy[:,:,-1]**2+dAdz[:,:,-1]**2)
                ynormal = dCdy[:,:,-1]/numpy.sqrt(dBdx[:,:,-1]**2+dBdy[:,:,-1]**2+dBdz[:,:,-1]**2)
                znormal = dCdz[:,:,-1]/numpy.sqrt(dCdx[:,:,-1]**2+dCdy[:,:,-1]**2+dCdz[:,:,-1]**2)

            self.BCdata[dataKey] = [xnormal,ynormal,znormal]/numpy.sqrt(xnormal**2+ynormal**2+znormal**2)


        return xnormal, ynormal, znormal



    def Reimann(self,direction,Rhoref,Uref,Vref,Wref,Pref,gamma,rhoi,Ui,Vi,Wi,Pi):
        """
        Riemann invariants
        index
        b - variables at the boundary index=0,ax,ay,az
        i - variables inside the domain index=1,n-2,...
        o - variables at outiside of the domain (inf)
        """

        # itializing free stream vectors
        rhoo = numpy.zeros((Ui.shape[0],1)) + Rhoref
        Uo   = numpy.zeros((Ui.shape[0],1)) + Uref
        Vo   = numpy.zeros((Ui.shape[0],1)) + Vref
        Wo   = numpy.zeros((Ui.shape[0],1)) + Wref
        Po   = numpy.zeros((Ui.shape[0],1)) + Pref

        # getting normal vectors
        nx,ny,nz = self.getnormal(direction)

        Vni = Ui*nx+Vi*ny+Wi*nz
        Vno = Uo*nx+Vo*ny+Wo*nz
        SSo = numpy.sqrt(gamma*Po/rhoo)   
        SSi = numpy.sqrt(gamma*Pi/rhoi)
        tmp = numpy.sqrt(Ui**2+Vi**2+Wi**2)/SSi   

        # Compute Riemann invariants for supersonic and subsonic cases
        Rplus  = numpy.where( tmp >= 1.0 ,  Vno + 2.0*SSo/(gamma -1.0) , Vni + 2.0*SSi/(gamma -1.0) )
        Rminus = numpy.where( tmp >= 1.0 ,  Vni - 2.0*SSi/(gamma -1.0) , Vno - 2.0*SSo/(gamma -1.0) )
        
        Vnormal = (Rminus + Rplus)*0.5
        SSb = (Rplus - Rminus)*(gamma-1.0)*0.25
        Ub = numpy.where( Vnormal >= 0.0 ,  Ui+(Vnormal-Vni)*nx , Uo+(Vnormal-Vno)*nx )
        Vb = numpy.where( Vnormal >= 0.0 ,  Vi+(Vnormal-Vni)*ny , Vo+(Vnormal-Vno)*ny )
        Wb = numpy.where( Vnormal >= 0.0 ,  Wi+(Vnormal-Vni)*nz , Wo+(Vnormal-Vno)*nz )

        aux = 1.0/(gamma-1.0)
        rhob = numpy.where( Vnormal >= 0.0 , ((rhoi**gamma*SSb**2)/(gamma*Pi))**aux , ((rhoo**gamma*SSb**2)/(gamma*Po))**aux)

        Pb  = rhob*SSb**2/gamma
        Etb  = Pb/(gamma -1.0) + 0.5*rhob*(Ub*Ub+Vb*Vb+Wb*Wb)

        return rhob,Ub,Vb,Wb,Etb,Pb




    def farfieldbc(self,direction):
        
        if type(direction) != type([]):
            direction = [direction]
        
        for d in direction:
            self.farfield( d )


    def farfield(self,direction):

        refdata = self.BCdata['farfield-properties-%s' % direction]
        Rhoref = refdata['rho0']
        Uref = refdata['u0']
        Vref = refdata['v0']
        Wref = refdata['w0']
        Pref = refdata['p0']
        gamma = refdata['gamma']

        rho = refdata['rho']
        u = refdata['u']
        v = refdata['v']
        w = refdata['w']
        p = refdata['p']


        # Direction switch
        if direction == 'x1':
            if self.pyranda.PyMPI.x1proc:
                rhoi = self.pyranda.variables[rho].data[1,:,:]
                Ui   = self.pyranda.variables[u].data[1,:,:]           
                Vi   = self.pyranda.variables[v].data[1,:,:]
                Wi   = self.pyranda.variables[w].data[1,:,:]
                Pi   = self.pyranda.variables[p].data[1,:,:]

                [rhob,Ub,Vb,Wb,Etb,Pb]=self.Reimann(direction,Rhoref,Uref,Vref,Wref,Pref,gamma,rhoi,Ui,Vi,Wi,Pi)

                # update variables
                self.pyranda.variables[rho].data[0,:,:] = rhob
                self.pyranda.variables[u].data[0,:,:] = Ub
                self.pyranda.variables[v].data[0,:,:] = Vb
                self.pyranda.variables[w].data[0,:,:] = Wb
                self.pyranda.variables[p].data[0,:,:] = Pb

                

        if direction == 'xn':
            if self.pyranda.PyMPI.xnproc:
                rhoi = self.pyranda.variables[rho].data[-2,:,:]
                Ui   = self.pyranda.variables[u].data[-2,:,:]           
                Vi   = self.pyranda.variables[v].data[-2,:,:]
                Wi   = self.pyranda.variables[w].data[-2,:,:]
                Pi   = self.pyranda.variables[p].data[-2,:,:]

                [rhob,Ub,Vb,Wb,Etb,Pb]=self.Reimann(direction,Rhoref,Uref,Vref,Wref,Pref,gamma,rhoi,Ui,Vi,Wi,Pi)

                # update variables
                self.pyranda.variables[rho].data[-1,:,:] = rhob
                self.pyranda.variables[u].data[-1,:,:] = Ub
                self.pyranda.variables[v].data[-1,:,:] = Vb
                self.pyranda.variables[w].data[-1,:,:] = Wb
                self.pyranda.variables[p].data[-1,:,:] = Pb


        if direction == 'y1':
            if self.pyranda.PyMPI.y1proc:
                rhoi = self.pyranda.variables[rho].data[:,1,:]
                Ui   = self.pyranda.variables[u].data[:,1,:]           
                Vi   = self.pyranda.variables[v].data[:,1,:]
                Wi   = self.pyranda.variables[w].data[:,1,:]
                Pi   = self.pyranda.variables[p].data[:,1,:]

                [rhob,Ub,Vb,Wb,Etb,Pb]=self.Reimann(direction,Rhoref,Uref,Vref,Wref,Pref,gamma,rhoi,Ui,Vi,Wi,Pi)

                # update variables
                self.pyranda.variables[rho].data[:,0,:] = rhob
                self.pyranda.variables[u].data[:,0,:] = Ub
                self.pyranda.variables[v].data[:,0,:] = Vb
                self.pyranda.variables[w].data[:,0,:] = Wb
                self.pyranda.variables[p].data[:,0,:] = Pb
                

        if direction == 'yn':
            if self.pyranda.PyMPI.ynproc:
                rhoi = self.pyranda.variables[rho].data[:,-2,:]
                Ui   = self.pyranda.variables[u].data[:,-2,:]
                Vi   = self.pyranda.variables[v].data[:,-2,:]
                Wi   = self.pyranda.variables[w].data[:,-2,:]
                Pi   = self.pyranda.variables[p].data[:,-2,:]

                [rhob,Ub,Vb,Wb,Etb,Pb]=self.Reimann(direction,Rhoref,Uref,Vref,Wref,Pref,gamma,rhoi,Ui,Vi,Wi,Pi)

                # update variables
                self.pyranda.variables[rho].data[:,-1,:] = rhob
                self.pyranda.variables[u].data[:,-1,:] = Ub
                self.pyranda.variables[v].data[:,-1,:] = Vb
                self.pyranda.variables[w].data[:,-1,:] = Wb
                self.pyranda.variables[p].data[:,-1,:] = Pb                                 


        if direction == 'z1':
            if self.pyranda.PyMPI.z1proc:
                rhoi = self.pyranda.variables[rho].data[:,:,1]
                Ui   = self.pyranda.variables[u].data[:,:,1]           
                Vi   = self.pyranda.variables[v].data[:,:,1]
                Wi   = self.pyranda.variables[w].data[:,:,1]
                Pi   = self.pyranda.variables[p].data[:,:,1]

                [rhob,Ub,Vb,Wb,Etb,Pb]=self.Reimann(direction,Rhoref,Uref,Vref,Wref,Pref,gamma,rhoi,Ui,Vi,Wi,Pi)

                # update variables
                self.pyranda.variables[rho].data[:,:,0] = rhob
                self.pyranda.variables[u].data[:,:,0] = Ub
                self.pyranda.variables[v].data[:,:,0] = Vb
                self.pyranda.variables[w].data[:,:,0] = Wb
                self.pyranda.variables[p].data[:,:,0] = Pb                  
                

        if direction == 'zn':
            if self.pyranda.PyMPI.znproc:
                rhoi = self.pyranda.variables[rho].data[:,:,-2]
                Ui   = self.pyranda.variables[u].data[:,:,-2]           
                Vi   = self.pyranda.variables[v].data[:,:,-2]
                Wi   = self.pyranda.variables[w].data[:,:,-2]
                Pi   = self.pyranda.variables[p].data[:,:,-2] 

                [rhob,Ub,Vb,Wb,Etb,Pb]=self.Reimann(direction,Rhoref,Uref,Vref,Wref,Pref,gamma,rhoi,Ui,Vi,Wi,Pi)

                # update variables
                self.pyranda.variables[rho].data[:,:,-1] = rhob
                self.pyranda.variables[u].data[:,:,-1] = Ub
                self.pyranda.variables[v].data[:,:,-1] = Vb
                self.pyranda.variables[w].data[:,:,-1] = Wb
                self.pyranda.variables[p].data[:,:,-1] = Pb    
