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
from __future__ import print_function
from mpi4py import MPI
import numpy 
import re
import sys,os
import time
import matplotlib.pyplot as plt

from .pyrandaMPI   import pyrandaMPI
from .pyrandaVar   import pyrandaVar
from .pyrandaEq    import pyrandaEq
from .pyrandaMesh  import pyrandaMesh
from .pyrandaIO    import pyrandaIO
from .pyrandaUtils import *
from .pyrandaTex   import pyrandaTex

                                              
class pyrandaSim:

    def __init__(self,name,meshOptions,silent=False):

        self.name = name
        self.silent = silent
        self.mesh = pyrandaMesh()

        # Parse shorthand mesh description if string
        if type(meshOptions) == type(''):
            self.mesh.makeMeshStr( meshOptions )
            meshOptions = self.mesh.options

        self.meshOptions = meshOptions

        if 'coordsys' in meshOptions:
            self.mesh.coordsys = meshOptions['coordsys']
        else:            
            raise ValueError('No suitable mesh type specified.')

        nx = meshOptions['nn'][0]
        ny = meshOptions['nn'][1]
        nz = meshOptions['nn'][2]

        dx = (meshOptions['xn'][0]-meshOptions['x1'][0])/max(nx-1,1)
        dy = (meshOptions['xn'][1]-meshOptions['x1'][1])/max(ny-1,1)
        dz = (meshOptions['xn'][2]-meshOptions['x1'][2])/max(nz-1,1)
        periodic = meshOptions['periodic']                    

        self.dx = dx
        self.dy = dy
        self.dz = dz
                
        self.nx = nx
        self.ny = ny
        self.nz = nz
        
        self.mesh.options = meshOptions
        self.mesh.dims  = meshOptions['dim']

        self.PyMPI = pyrandaMPI( self.mesh )
        self.mesh.PyMPI = self.PyMPI

        self.mesh.makeMesh()

        # IO setup
        self.PyIO = pyrandaIO( self.name, self.PyMPI )
        
        
        # Grab some scalars off of parcop.mesh
        self.zero = self.PyMPI.emptyScalar()
        
        
        self.equations = []
        self.conserved = []
        self.variables = {}
        self.initial_conditions = []
        self.nPDE = 0
        self.nALG = 0

        # Package info
        self.packages = {}

        # Compute sMap
        self.get_sMap()

        # Time
        self.time = 0.0
        self.cycle = 0
        self.vizDumpHistory = []

        
        # Print startup message
        self.iprint( code() )
        self.iprint( version() )
        self.iprint( icopyright() )
        

    def iprint(self,sprnt):
        if self.PyMPI.master and (not self.silent):
            print(sprnt)
            sys.stdout.flush() 
        
    def addPackage(self,package):
        
        self.packages[package.name] = package
        package.get_sMap()
        self.sMap.update( package.sMap )
            
    def allocate(self):
        """
        Loop over vars and allocate the data
        """
        for v in self.variables:
            self.variables[v].__allocate__(self.PyMPI)
        
    
    def var(self,name):
        if name in self.variables.keys():
            return self.variables[name]
        else:
            raise ValueError('Error: variable name: %s not found in database' % name)


    def addVar(self,name,kind=None,rank='scalar'):
        if name not in self.variables:
            var = pyrandaVar(name,kind,rank)
            self.variables[name] = var

    def addEqu(self,equation):

        peq = pyrandaEq(equation,self.sMap,self.PyMPI) 
        self.equations.append( peq )
        if peq.kind == 'PDE':
            self.nPDE += 1
            self.conserved.append( peq.LHS[0] )
            if len( peq.LHS ) > 1:
                print('Warning... only single return values for PDEs allowed')
                exit()
        elif peq.kind == 'ALG':
            self.nALG += 1
        else:
            raise ValueError('No suitable equation type specified for %s' % peq.eqstr)

    def EOM(self,eom):
        """
        Higher level wrapper to make equations of motion from a single string
        Will add eqautions and variables as needed.
        """

        # Split up the equation lines
        eom_lines = filter(None,eom.split('\n'))
        eom_lines = [el for el in eom_lines if el.strip()[0] != '#']  # Comments work
        var_names = []

        #### Get unique set of variables - scalars ####
        for eq in eom_lines:
            evars = findVar( eq,'scalar')
            var_names += evars
        var_names = list(set(var_names))
        for evar in var_names:
            self.addVar(evar,kind='conserved')   # Todo: classify variables

        #### Get unique set of variables - vectors ####
        var_names = []
        for eq in eom_lines:
            evars = findVar( eq,'vector')
            var_names += evars
        var_names = list(set(var_names))
        for evar in var_names:
            self.addVar(evar,rank='vector',kind='conserved')   
                                
        for variable in self.variables:
            self.iprint('Adding variables: %s' %  variable )

        #### Add equations of motion ####
        self.iprint('Adding equations of motion: ')

        eqN = 1
        for eq in eom_lines:
            self.addEqu( eq )
            self.iprint( '(%s)   %s' % (eqN,eq) )
            eqN += 1
            

        # Set up the variables in memory
        self.allocate()
            
    def checkForNan(self,names=[]):

        nans = False
        svars = ''
        if not names:
            names = self.variables
        for ivar in names:
            try:
                myvar = self.variables[ivar]
                if numpy.isnan(myvar.data).any():
                    nans = True
                    svars += ivar + ' '
            except:
                print("%s is not a variable" % ivar)
                #import pdb
                #pdb.set_trace()
        return svars

    def setIC(self,ics):
        """
        Evaluate the initial conditions and then update variables
        """
        # Split up the equation lines
        ic_lines = filter(None,ics.split('\n'))
        ic_lines = [el.replace(' ','') for el in ic_lines ]  # Comments work
        ic_lines = filter(None,ic_lines)
        ic_lines = [el for el in ic_lines if el.strip()[0] != '#']  # Comments work

        var_names = []

        #### Get unique set of variables ####
        for eq in ic_lines:
            evars = findVar( eq,'scalar')
            var_names += evars

        var_names = list(set(var_names))

        for evar in var_names:
            self.addVar(evar,kind='conserved')   # Todo: classify variables
        
        # Actually compute the Initial Conditions
        for ic in ic_lines:
            ic_mod = ic #+ '+self.emptyScalar()'
            exec(fortran3d(ic_mod,self.sMap))
        
        for ic in ic_lines:
            self.initial_conditions.append( ic )
            
        # Check for nans here
        snans = self.checkForNan( var_names )
        if snans:
            self.iprint("Found some nans in inits: %s" % snans)
            exit()

        self.updateVars()
        snans = self.checkForNan()
        if snans:
            self.iprint("Found some nans in Update after init: %s" % snans)
            exit()


        
    def emptyScalar(self):
        return self.PyMPI.emptyScalar()
        
    def updateFlux(self):
        #
        ncons = self.nPDE
        shape = self.mesh.shape[:]
        shape.append( ncons )        
        # Only ddt() terms are extracted
        flux = {} #numpy.asfortranarray( numpy.zeros( shape ) )
        #ieq = 0
        for eqo in self.equations:
            eq = eqo.eqstr
            if ( eqo.kind == 'PDE' ):
                lhs = eqo.LHS[0]
                Srhs = eq.split('=')[1]  # This is a string to evaluate
                flux[lhs] = eqo.RHS(self)

        return flux

    def updateVars(self):
        #
        # Update the equations
        #
        for eq in self.equations:

            if not eq.active:
                continue

            if ( eq.kind == 'ALG'):            
                rhs = eq.RHS(self)

                # If no LHS , then on to the next eq
                if not eq.LHS:
                    continue
                
                if eq.rank == 1:
                    self.variables[eq.LHS[0]].data = rhs
                else:
                    for ii in range(len(rhs)):
                        #import pdb
                        #pdb.set_trace()
                        self.variables[eq.LHS[ii]].data = rhs[ii]
       

        
        
    def write(self,wVars=[]):
        """ 
        Write viz file 
        """
        procInt = 5
        visInt  = 5
        if not wVars:
            wVars = self.conserved

        dumpFile = 'vis' + str(self.cycle).zfill(visInt)
        
        if self.PyMPI.master == 1:
            try:
                os.mkdir(os.path.join(self.PyIO.rootname, dumpFile))
            except:
                pass

        self.PyMPI.comm.barrier()   # Wait for directory to be made
        rank = self.PyMPI.comm.rank
        dumpFile = os.path.join( self.PyIO.rootname,
                                 dumpFile,
                                 'proc-%s.%s' % (str(rank).zfill(procInt),str(self.cycle).zfill(visInt)))


        suff = 'vtk'
        self.PyIO.makeDumpVTK(self.mesh,self.variables,wVars,dumpFile)

        self.vizDumpHistory.append( [self.cycle, self.time] )

        # Write .visit file
        if self.PyMPI.master == 1:
            vid = open( os.path.join(self.PyIO.rootname, 'pyranda.visit' ) , 'w')
            vid.write("!NBLOCKS %s \n" % self.PyMPI.comm.Get_size() )
            for vdump in self.vizDumpHistory:
                iv = vdump[0]
                for p in range(self.PyMPI.comm.Get_size()):
                    vid.write("%s\n" % os.path.join('vis' + str(iv).zfill(visInt),
                                                    'proc-%s.%s.%s' % (str(p).zfill(procInt),str(iv).zfill(visInt),suff ) ) )
            vid.close()
                    
            
        

        

    def writeGrid(self):
        """ 
        Write grid file 
        """
        wlen = 3
        shape = list(self.mesh.shape)
        shape.append( wlen )
        iodata = numpy.zeros( shape )
        for i in range( 3 ):
            iodata[:,:,:,i] = self.mesh.coords[i]

        dumpName = 'grid' 
        self.PyIO.makeDump(iodata,dumpName)

        
        

        
                        
    def ddx(self,val):
        if self.nx <= 1:
            return 0.0
        return self.PyMPI.der.ddx( val )

    def ddy(self,val):
        if self.ny <= 1:
            return 0.0
        return self.PyMPI.der.ddy( val )

    def ddz(self,val):
        if self.nz <= 1:
            return 0.0
        return self.PyMPI.der.ddz( val )

        
    def div(self,f1,f2=None,f3=None):

        if (type(f2) == type(None) and type(f3) == type(None)):
            if self.nx > 1:
                return self.PyMPI.der.div(f1,self.zero,self.zero)
            if self.ny > 1:
                return self.PyMPI.der.div(self.zero,f1,self.zero)
            if self.nz > 1:
                return self.PyMPI.der.div(self.zero,self.zero,f1)

        elif type(f3) == type(None):
            if (self.nx > 1) and (self.ny > 1):
                return self.PyMPI.der.div(f1,f2,self.zero)
            if (self.nz > 1) and (self.ny > 1):
                return self.PyMPI.der.div(self.zero,f1,f2)
            if (self.nz > 1) and (self.nx > 1):
                return self.PyMPI.der.div(f1,self.zero,f2)
        else:            
            return self.PyMPI.der.div(f1,f2,f3)
    
    def grad(self,val):
        return [self.ddx(val),self.ddy(val),self.ddz(val)]

    def laplacian(self,val):
        return self.PyMPI.der.laplacian( val )

    def ring(self,val):
        return self.PyMPI.der.ring( val )

    def filterx(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_x(val, f_tilde)
        return f_tilde

    def filtery(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_y(val, f_tilde)
        return f_tilde

    def filterz(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.fil.filter_z(val, f_tilde)
        return f_tilde 


    def filter(self,val):
        return self.PyMPI.fil.filter(val)

    def filterOld(self,val):
        if self.nx > 1:
            f1 = self.filterx(val)
        else:
            f1 = val
        if self.ny > 1:
            f2 = self.filtery(f1)
        else:
            f2 = f1
        if self.nz > 1:
            f1 = self.filterz(f2)
        else:
            f1 = f2
        return f1

    def gfilterx(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.gfil.filter_x(val, f_tilde)
        return f_tilde

    def gfiltery(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.gfil.filter_y(val, f_tilde)
        return f_tilde

    def gfilterz(self,val):
        f_tilde = self.emptyScalar()
        self.PyMPI.gfil.filter_z(val, f_tilde)
        return f_tilde 

    def gfilter(self,val):
        return self.PyMPI.gfil.filter(val)
    
    def gfilterOld(self,val):
        if self.nx > 1:
            f1 = self.gfilterx(val)
        else:
            f1 = val
        if self.ny > 1:            
            f2 = self.gfiltery(f1)
        else:
            f2 = f1
        if self.nz > 1:
            f1 = self.gfilterz(f2)
        else:
            f1 = f2
        return f1

    
    def rk4(self,time,dt):

        Ark = [0.0]*5
        Ark[0] = 0.0;
        Ark[1] = -6234157559845./12983515589748.;
        Ark[2] = -6194124222391./4410992767914.;
        Ark[3] = -31623096876824./15682348800105.;
        Ark[4] = -12251185447671./11596622555746.;

        Brk = [0.0]*5
        Brk[0] = 494393426753./4806282396855.;
        Brk[1] = 4047970641027./5463924506627.;
        Brk[2] = 9795748752853./13190207949281.;
        Brk[3] = 4009051133189./8539092990294.;
        Brk[4] = 1348533437543./7166442652324.;

        eta = [0.0]*5
        eta[0] = 494393426753./4806282396855.;
        eta[1] = 4702696611523./9636871101405.;
        eta[2] = 3614488396635./5249666457482.;
        eta[3] = 9766892798963./10823461281321.;
        eta[4] = 1.0;

        #	Initialize some intermediate arrays
        ncons = self.nPDE
        shape = tuple(self.mesh.shape)

        tmp1 = {}
        tmp2 = {}
        PHI  = {}
        for U in self.conserved: 
            tmp1[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            tmp2[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            PHI[U]  = numpy.asfortranarray( numpy.zeros(  shape ) )
        
        # Get primative flow variables
        #self.updateVars()
        time_i = time
        #import pdb
        #pdb.set_trace()
        for ii in range(5):
            #    ii
            FLUX = self.updateFlux()
            for U in self.conserved:
                tmp1[U] =  Ark[ii]*PHI[U]
                PHI[U]  =  dt*FLUX[U] + tmp1[U]
                tmp2[U] =  Brk[ii]*PHI[U]
                self.variables[U].data =  self.variables[U].data + tmp2[U]
            time = time_i + eta[ii]*dt
            self.time = time
            self.updateVars()

        self.cycle += 1

        return time

    def get_sMap(self):
        sMap = {}
        
        #sMap["div(#arg#)"] = ""
        #if self.nx > 1:
        #    sMap["div(#arg#)"] += "self.ddx(#arg#[:,:,:,0])"
        #if self.ny > 1:
        #    sMap["div(#arg#)"] += "+self.ddy(#arg#[:,:,:,1])"
        #if self.nz > 1:
        #    sMap["div(#arg#)"] += "+self.ddz(#arg#[:,:,:,2])"

        # Simple find/replace mappings
        sMap['div(' ] = 'self.div('
        sMap['ddx(' ] = 'self.ddx('
        sMap['ddy(' ] = 'self.ddy('
        sMap['ddz(' ] = 'self.ddz('
        sMap['fbar('] = 'self.filter('
        sMap['gbar('] = 'self.gfilter('
        sMap['grad('] = 'self.grad('
        sMap['simtime'] = 'self.time'
        sMap['lap(' ] = 'self.laplacian('
        sMap['ring(' ] = 'self.ring('
        sMap['sum(' ] = 'self.PyMPI.sum3D('
        sMap['sign(' ] = 'numpy.sign('
        sMap['dot(' ] = 'numpy.dot('
        sMap['abs(' ] = 'numpy.abs('
        sMap['sqrt(' ] = 'numpy.sqrt('
        sMap['sin('] = 'numpy.sin('
        sMap['cos('] = 'numpy.cos('
        sMap['tanh('] = 'numpy.tanh('
        sMap['exp('] = 'numpy.exp(' 
        sMap['where('] = 'numpy.where('
        sMap['max('] = 'numpy.maximum('
        sMap['min('] = 'numpy.minimum('
        sMap['3d()'] = 'self.emptyScalar()'
        sMap[':pi:'] = 'numpy.pi'
        
        sMap['meshx']   = 'self.mesh.coords[0]'
        sMap['meshy']   = 'self.mesh.coords[1]'
        sMap['meshz']   = 'self.mesh.coords[2]'
        self.sMap = sMap
        
    def euler(self,time,dt):

        #	Initialize some intermediate arrays
        ncons = self.nPDE
        shape = tuple(self.mesh.shape)

        # Get primative flow variables
        time_i = time
        FLUX = self.updateFlux()
        for U in self.conserved:
            self.variables[U].data += FLUX[U]*dt  

        time = time_i + dt
        self.time = time
        self.updateVars()
        return time

    def rk4_step(self,dt,ii,PHI):

        Ark = [0.0]*5
        Ark[0] = 0.0;
        Ark[1] = -6234157559845./12983515589748.;
        Ark[2] = -6194124222391./4410992767914.;
        Ark[3] = -31623096876824./15682348800105.;
        Ark[4] = -12251185447671./11596622555746.;

        Brk = [0.0]*5
        Brk[0] = 494393426753./4806282396855.;
        Brk[1] = 4047970641027./5463924506627.;
        Brk[2] = 9795748752853./13190207949281.;
        Brk[3] = 4009051133189./8539092990294.;
        Brk[4] = 1348533437543./7166442652324.;

        eta = [0.0]*5
        eta[0] = 494393426753./4806282396855.;
        eta[1] = 4702696611523./9636871101405.;
        eta[2] = 3614488396635./5249666457482.;
        eta[3] = 9766892798963./10823461281321.;
        eta[4] = 1.0;

        #	Initialize some intermediate arrays
        ncons = self.nPDE
        shape = tuple(self.mesh.shape)

        tmp1 = {}
        tmp2 = {}

        for U in self.conserved: 
            tmp1[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            tmp2[U] = numpy.asfortranarray( numpy.zeros(  shape ) )
            if not PHI:
                PHI  = {}
                PHI[U]  = numpy.asfortranarray( numpy.zeros(  shape ) )
        
        # Get primative flow variables

        #for ii in range(5):
            #    ii
        FLUX = self.updateFlux()
        for U in self.conserved:
            tmp1[U] =  Ark[ii]*PHI[U]
            PHI[U]  =  dt*FLUX[U] + tmp1[U]
            tmp2[U] =  Brk[ii]*PHI[U]
            self.variables[U].data += tmp2[U]
        self.updateVars()        
        return PHI

    

    def get_tMap(self):
        tMap = {}
        
        # Simple find/replace mappings
        tMap['ddx(#arg#)' ] = r'\frac{\partial #arg#}{\partial x}'
        tMap['ddy(#arg#)' ] = r'\frac{\partial #arg#}{\partial y}'
        tMap['ddz(#arg#)' ] = r'\frac{\partial #arg#}{\partial z}'
        tMap['ddt(#arg#)' ] = r'\frac{\partial #arg#}{\partial t}'

        tMap['meshx'] = 'x'
        tMap['meshy'] = 'y'
        tMap['meshz'] = 'z'

        tMap[':pi:'] = r'\pi '

        #tMap['exp(#arg#)'] = r'\exp( #arg# )'
        #tMap['sqrt(#arg#)'] = r'\sqrt( #arg# )'
        tMap['exp('] = r'\exp('
        tMap['sqrt('] = r'\sqrt('

        tMap['*' ] = ''
        tMap['='] = r'&='
        

        return tMap


    def setupLatex(self):

        self.latex = pyrandaTex(self)
        
