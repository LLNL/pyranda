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
import os
import numpy 



class pyrandaIO:
    """
    Class to handle IO of pyranda objects and files
    """
    
    def __init__(self,rootname,pympi):
        
        self.rootname = rootname
        self.PyMPI = pympi
        
        self.variables = []  # List of variables to be included
        
        if self.PyMPI.master == 1:
            try:
                os.mkdir(rootname)
            except:
                pass
        self.PyMPI.comm.barrier()   # Wait for directory to be made
            
    def makeDump(self,data,dumpName):
        
        dumpFile = dumpName
        if self.PyMPI.master == 1:
            try:
                os.mkdir(os.path.join(self.rootname, dumpFile))
            except:
                pass

        rank = self.PyMPI.comm.rank
        file_name = os.path.join( self.rootname, dumpFile, 'p' + str(rank).zfill(8) )


        fd = open(file_name, 'wb')            
        #fwrite(fd, data.size, data)
        data.tofile(fd)
        fd.close()


    def makeDumpVTK(self,mesh,variables,varList,dumpFile):

        ghost = True

        if ghost:
            fx = self.PyMPI.ghost(mesh.coords[0].data[:,:,:] )
            fy = self.PyMPI.ghost(mesh.coords[1].data[:,:,:] )
            fz = self.PyMPI.ghost(mesh.coords[2].data[:,:,:] )
        else:
            fx = mesh.coords[0].data[:,:,:]
            fy = mesh.coords[1].data[:,:,:]
            fz = mesh.coords[2].data[:,:,:]
            
        ax = fx.shape[0]
        ay = fx.shape[1]
        az = fx.shape[2]
        
        fx = fx.flatten(order='F')
        fy = fy.flatten(order='F')
        fz = fz.flatten(order='F')
        
        mesh = numpy.vstack((fx,fy,fz)).transpose()
        
        fid = open(dumpFile + '.vtk','w')
                
        fid.write("# vtk DataFile Version 3.0 \n")
        fid.write("vtk output \n")
        fid.write("ASCII \n")
        fid.write("DATASET STRUCTURED_GRID \n")
        fid.write("DIMENSIONS  %s %s %s  \n" % (ax,ay,az))
        fid.write("POINTS %s float  \n" % len(fx))

        numpy.savetxt(fid,mesh,fmt='%f')

        fid.write("POINT_DATA %s  \n" % len(fx) )

        for var in varList:
            fid.write("SCALARS %s float \n" % var)
            fid.write("LOOKUP_TABLE default \n")
            if ghost:
                gdata = self.PyMPI.ghost( variables[var].data )
            else:
                gdata = variables[var].data
            numpy.savetxt(fid,gdata.flatten(order='F') ,fmt='%f')
        
        fid.close()
        

    def makeDumpTec(self,mesh,variables,varList,dumpFile):
        
        fx = mesh.coords[0].data[:,:,:].flatten()
        fy = mesh.coords[1].data[:,:,:].flatten()
        fz = mesh.coords[2].data[:,:,:].flatten()
        
        fid = open(dumpFile + '.tec','w')
        
        ax = self.PyMPI.chunk_3d_size[0]
        ay = self.PyMPI.chunk_3d_size[1]
        az = self.PyMPI.chunk_3d_size[2]

        fid.write('VARIABLES = "X", "Y", "Z" ')
        for var in varList:
            fid.write(', "%s" ' % var)
        fid.write('\n')
                
        fid.write("ZONE  I=%s, J=%s, K=%s, DATAPACKING=BLOCK  \n" % (ax,ay,az))

        # Mesh
        numpy.savetxt(fid,fx,fmt='%f')
        numpy.savetxt(fid,fy,fmt='%f')
        numpy.savetxt(fid,fz,fmt='%f')

        # Variables
        for var in varList:
            numpy.savetxt(fid,variables[var].data.flatten() ,fmt='%f')
        
        fid.close()

        
        
    def makeMIR(self):
        

        # Write meta file for viz
        form ="""
VERSION 2.0
zonal: yes                # Zonal or Nodal?
curvilinear: yes          # Read grid files?
gridfiles: grid/p%08d     # grid files
datafiles: vis%07d/p%08d  # data files
fileorder: XYZ # processor order
domainsize:       #AX#    #AY#    #AZ# # overall nodal dimensions
blocksize:        #AX#    #AY#    #AZ# # nodal dimensions per miranda processor
interiorsize:     #AX#    #AY#    #AZ# # nodal dimensions per interior processor
bndrysize:        #AX#    #AY#    #AZ# # nodal dimensions per boundary processor
origin:    #X1#  #Y1#  #Z1#            # xmin, ymin, zmin for rectilinear
spacing:   #DX#  #DY#  #DZ#            # dx, dy, dz for rectilinear
variables:   #NVARS#  # number of variables
  #VARS#
timesteps:    #NVIZ#  # number of times to plot
  #CYCTIME#
"""     
        #form = form.replace('#AX#',str(self.PyMPI.chunk_3d_size[0]))
        #form = form.replace('#AY#',str(self.PyMPI.chunk_3d_size[1]))
        #form = form.replace('#AZ#',str(self.PyMPI.chunk_3d_size[2]))

        #form = form.replace('#X1#',str(self.mesh.options['x1'][0]))
        #form = form.replace('#Y1#',str(self.mesh.options['x1'][1]))
        #form = form.replace('#Z1#',str(self.mesh.options['x1'][2]))
        
        #form = form.replace('#DX#',str(self.dx))
        #form = form.replace('#DY#',str(self.dy))
        #form = form.replace('#DZ#',str(self.dz))

        #form = form.replace('#NVARS#',str(wlen))
        #svars = ''
        #for ivar in wVars:
        #    svars += "  %s 1 %f %f \n" % (ivar,
        #                                  self.PyMPI.min3D( self.variables[ivar].data),
        #                                  self.PyMPI.max3D( self.variables[ivar].data) )
        #form = form.replace("#VARS#",svars)
        
        # Time history
        #form = form.replace("#NVIZ#",str(len( self.vizDumpHistory) ))
        #sviz = ''
        #for iv in self.vizDumpHistory:
        #    sviz += "  %s  %f \n" % ( str(iv[0]).zfill(7), iv[1] )
        #form = form.replace("#CYCTIME#",sviz)

        #if self.PyMPI.master:
        #    mfid = open( os.path.join( self.PyIO.rootname, 'pyranda.mir'),'w+')
        #    mfid.write(form)
        #    mfid.close()
        


