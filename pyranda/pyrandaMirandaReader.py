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
import sys,os
from os.path import join
from .pyranda import pyrandaSim

# Gaurd against hdf5 dependency
try:
    import h5py
except:
    print("h5py is required for mirandaReader to work")
    pass

# Static paths for sidre/blueprint import
path_proc  = '/sidre/groups/marbl/groups/miranda/groups/blueprint/views/'
path_coord = '/sidre/groups/marbl/groups/miranda/groups/blueprint/groups/coordsets/groups/coords/groups/dims/views/'
path_grid  = '/sidre/groups/marbl/groups/miranda/groups/blueprint/groups/coordsets/groups/coords/groups/spacing/views/'
path_field = '/sidre/external/marbl/miranda/blueprint/fields/'
path_mats  = '/sidre/external/marbl/miranda/blueprint/matsets/mats/volume_fractions/'

class mirandaReader:

    def __init__(self, rootdir, var_list, periodic=[False,False,False]):
        
        self.rootdir = rootdir
        self.var_list = var_list
        
        # Rootfile is directory where the marbl_* files are located
        cycle = 0
        self.data_dir = join(rootdir,'marbl_' + str(cycle).zfill(7))

        # Check if file exists
        if not os.path.isdir( self.data_dir ):
            print("Error: Cant read restart file %s" % self.data_dir )
            return None
        
        # Read in domain meta-data from master proc
        proc_num = 0
        filename = join(self.data_dir,'marbl_' + str(cycle).zfill(7) + '_' + str(proc_num).zfill(7) + '.hdf5')
        header = 'datagroup_' + str(proc_num).zfill(7)

        # Load hdf5 file
        f = h5py.File(filename,'r')

        # Processor arrangement
        px = int(f[header + path_proc + 'px/value'].value)
        py = int(f[header + path_proc + 'py/value'].value)
        pz = int(f[header + path_proc + 'pz/value'].value)
        self.procs = [px,py,pz]

        # Grid spacing
        dx = f[header + path_grid + 'dx/value'].value
        dy = f[header + path_grid + 'dy/value'].value
        dz = f[header + path_grid + 'dz/value'].value

        # Import data as numpy array
        # Grab grid size for verificiation  - Static location in blueprint?
        ax = int(f[header + path_coord + 'i/value'].value) - 1
        ay = int(f[header + path_coord + 'j/value'].value) - 1
        az = int(f[header + path_coord + 'k/value'].value) - 1

        # Need to get a mesh object here
        x1,xn = 0, dx*ax*px
        y1,yn = 0, dy*ay*py
        z1,zn = 0, dz*az*pz
        nx,ny,nz = ax*px,ay*py,az*pz
        mesh_options = {}
        mesh_options['x1'] = [x1,y1,z1]
        mesh_options['xn'] = [xn,yn,zn]
        mesh_options['nn'] = [nx,ny,nz]
        mesh_options['periodic'] = periodic
        
        # Make pyranda instance
        self.pysim = pyrandaSim( rootdir + '_pyranda', mesh_options )

        # Add variables to pyranda-state
        for v in self.var_list:
            self.pysim.addVar(v)
        self.pysim.allocate()

        # Make a procmap for the source data/MPI
        self.procMap = {}
        ipp = 0
        for ipz in range(pz):
            for ipy in range(py):
                for ipx in range(px):
                        
                    gx1 = ax*ipx + 1 - 1
                    gy1 = ay*ipy + 1 - 1
                    gz1 = az*ipz + 1 - 1                          

                    gxn = ax*(ipx+1) - 1
                    gyn = ay*(ipy+1) - 1
                    gzn = az*(ipz+1) - 1
                    
                    self.procMap['%s-g1' % ipp] = numpy.array([gx1,gy1,gz1])
                    self.procMap['%s-gn' % ipp] = numpy.array([gxn,gyn,gzn])
                    
                    ipp += 1
                    
    def update(self, cycle):
        """
        Import Miranda data from specified cycle
        """
        
        # Update data directory to current cycle
        self.data_dir = join(self.rootdir,'marbl_' + str(cycle).zfill(7))

        # Check if file exists
        if not os.path.isdir( self.data_dir ):
            self.pysim.iprint("Error: Cant read restart file %s" % self.data_dir )
            return None

        readChunkMiranda( self.pysim, self.procs, self.procMap, self.data_dir, cycle, self.var_list)

        return self.pysim
        
def readChunkMiranda(pysim,procs,procMap,dump,cycle,var_list):
    """
    Same as readData but only reads in global range of data given by
    irange.
    """

    Rx = [0]*2
    Ry = [0]*2
    Rz = [0]*2

    # This procs extents
    Rx[0] = pysim.PyMPI.chunk_3d_lo[0]    #irange[0]
    Rx[1] = pysim.PyMPI.chunk_3d_hi[0]+1  #irange[1]
    Ry[0] = pysim.PyMPI.chunk_3d_lo[1]    #irange[2]
    Ry[1] = pysim.PyMPI.chunk_3d_hi[1]+1  #irange[3]
    Rz[0] = pysim.PyMPI.chunk_3d_lo[2]    #irange[4]
    Rz[1] = pysim.PyMPI.chunk_3d_hi[2]+1  #irange[5]

    # Restart deomain info
    nprocs = procs[0]*procs[1]*procs[2]
    ax = pysim.nx / procs[0]
    ay = pysim.ny / procs[1]
    az = pysim.nz / procs[2]
    nshape = (ax,ay,az,len(pysim.variables))
    

    for iproc in range(nprocs):

        g1 = procMap['%s-g1' % iproc] 
        gn = procMap['%s-gn' % iproc] + 1

        # Shift left point if node data
        iff = 0;jff = 0;kff = 0;

        c1 = (Rx[1] in range(g1[0],gn[0]) )
        c2 = (Rx[0] in range(g1[0],gn[0]) )
        c3 = ( (g1[0] and gn[0]) in range(Rx[0],Rx[1]+1) )
        CX = c1 or c2 or c3

        c1 = (Ry[1] in range(g1[1],gn[1]) )
        c2 = (Ry[0] in range(g1[1],gn[1]) )
        c3 = ( (g1[1] and gn[1]) in range(Ry[0],Ry[1]+1) )
        CY = c1 or c2 or c3

        c1 = (Rz[1] in range(g1[2],gn[2]) )
        c2 = (Rz[0] in range(g1[2],gn[2]) )
        c3 = ( (g1[2] and gn[2]) in range(Rz[0],Rz[1]+1) )
        CZ = c1 or c2 or c3

        if ( CX and CY and CZ ):

            Li1 = numpy.max( (0 , Rx[0] - g1[0] ) ) + iff
            Lif = numpy.min( (Rx[1] , gn[0] ) ) - g1[0] + iff
            Ki1 = numpy.max( (Rx[0] , g1[0]) ) - Rx[0]
            Kif = Ki1 + (Lif-Li1)

            Lj1 = numpy.max( (0 , Ry[0] - g1[1] ) ) + jff
            Ljf = numpy.min( (Ry[1] , gn[1] ) ) - g1[1] + jff
            Kj1 = numpy.max( (Ry[0] , g1[1]) ) - Ry[0]
            Kjf = Kj1 + (Ljf-Lj1)

            Lk1 = numpy.max( (0 , Rz[0] - g1[2] ) ) + kff
            Lkf = numpy.min( (Rz[1] , gn[2] ) ) - g1[2] + kff
            Kk1 = numpy.max( (Rz[0] , g1[2]) ) - Rz[0]
            Kkf = Kk1 + (Lkf-Lk1)

            # Read proc file
            proc_num = iproc          
            filename = join(dump,'marbl_' + 
                            str(cycle).zfill(7) + '_' + 
                            str(proc_num).zfill(7) + '.hdf5')

            # Load hdf5 file
            f = h5py.File(filename,'r')

            for var in var_list:
                DATA = import_field(f,iproc,var)    
                pysim.variables[var].data[Ki1:Kif,
                                          Kj1:Kjf,
                                          Kk1:Kkf] = DATA[Li1:Lif,
                                                          Lj1:Ljf,
                                                          Lk1:Lkf]
            
def convert(str_array):     # Convert numpy array of binary encoded chars to string
    s = ""
    for ele in str_array:
        s += ele.decode('UTF-8')
    return s

def import_field(data, proc_num, var):
    
    header = 'datagroup_' + str(proc_num).zfill(7)
    nx = int(data[header + path_coord + 'i/value'].value) - 1
    ny = int(data[header + path_coord + 'j/value'].value) - 1
    nz = int(data[header + path_coord + 'k/value'].value) - 1
    n = (nx,ny,nz)

    if var[0] == 'Y':
        tmp = data[header + path_mats + str(var[1:]).zfill(3)]       # Temp data holder
        arr = numpy.zeros(tmp.shape,dtype='float64')    # Initialize array using dataset shape
        tmp.read_direct(arr)       # Directly read in hdf5 dataset to numpy array      
        #print('Finished importing ' + var + str(spec_num).zfill(3) + '!')      
        
    else:
        tmp = data[header + path_field + var + '/values']       # Temp data holder
        arr = numpy.zeros(tmp.shape,dtype='float64')    # Initialize array using dataset shape
        tmp.read_direct(arr)       # Directly read in hdf5 dataset to numpy array
        #print('Finished importing ' + var + '!')
    
    return numpy.reshape(arr,n,order='F')  
