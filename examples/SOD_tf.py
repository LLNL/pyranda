import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC
from pyranda.pyranda import pyrandaRestart


#teach = True
teach = False


if not teach:
    # TENSOR FLOW
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers



## Define a mesh
L = numpy.pi * 2.0
Npts = 200
Lp = L * (Npts-1.0) / Npts

imesh = """
xdom = (0.0, Lp, Npts)
""".replace('Lp',str(Lp)).replace('Npts',str(Npts))

# Initialize a simulation object on a mesh
ss = pyrandaSim('sod',imesh)
ss.addPackage( pyrandaBC(ss) )


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) 
:beta:      =  gbar( ring(:div:) * :rho:) * 7.0e-2
:tau:       =  :beta:*:div:
:beta:      =  :mlbeta: * 1.0
:mlbeta:    =  :beta: * 0.0
:tmp:       =  :beta: * 0.0
# Apply constant BCs
bc.extrap(['rho','Et'],['x1'])
bc.const(['u'],['x1','xn'],0.0)
"""

# Add the EOM to the solver
ss.EOM(eom)


# Initial conditions SOD shock tube in 1d
ic = """
:gamma: = 1.4
:Et:  = gbar( where( meshx < pi, 1.0/(:gamma:-1.0) , .1 /(:gamma:-1.0) ) )
:rho: = gbar( where( meshx < pi, 1.0    , .125 ) )
"""

# Set the initial conditions
ss.setIC(ic)
    

# Write a time loop
time = 0.0

# Approx a max dt and stopping time
v = 1.0
dt_max = v / ss.mesh.nn[0] * 0.75
tt = L/v * .25 #dt_max

# Start time loop
dt = dt_max
cnt = 1
viz_freq = 25
pvar = 'beta'
viz = True


ML_data = numpy.zeros( (8,400 * (Npts-6)) )
mlcnt = 0

if not teach:
    new_model = tf.keras.models.load_model('my_first_model.h5')

while tt > time:

    # Update the EOM and get next dt
    #time = ss.rk4(time,dt)
    time = ss.euler(time,dt*.2)
    dt = min(dt_max, (tt - time) )

    cnt += 1
    
    # Compute ML data for bulk viscosity
    ml_beta = ss.var('mlbeta').data
    p_data = numpy.zeros( (Npts-6,7) )
    for i in range(3,Npts-3):
        if mlcnt < (400*(Npts-6)-1):
            
            ML_data[0,mlcnt] = ss.var('beta').data[i,0,0]
            ML_data[1,mlcnt] = ss.var('u').data[i-3,0,0]
            ML_data[2,mlcnt] = ss.var('u').data[i-2,0,0]
            ML_data[3,mlcnt] = ss.var('u').data[i-1,0,0]
            ML_data[4,mlcnt] = ss.var('u').data[i,0,0]
            ML_data[5,mlcnt] = ss.var('u').data[i+1,0,0]
            ML_data[6,mlcnt] = ss.var('u').data[i+2,0,0]
            ML_data[7,mlcnt] = ss.var('u').data[i+3,0,0]

            # Only teach near shocks
            div = ss.var('div').data[i-3:i+4,0,0]
            if min( div ) < -.1:                
                mlcnt += 1
                ss.var('tmp').data[i,0,0] = 1.0
            # Apply model to data
            if not teach:
                #if (cnt%viz_freq == 0):                
                u_mean = .271586
                u_std = .3912
                p_data[i-3,:] = (( ML_data[1:,mlcnt-1] - u_mean ) / u_std)

    if not teach:
        #if (cnt%viz_freq == 0):
        ml_beta[3:-3,0,0] = new_model.predict( p_data )[:,0]
            
        ss.var('mlbeta').data = ml_beta * ss.var('tmp').data

        # Print some output
    ss.iprint("%s -- %s" % (cnt,time)  )

    if viz:

        if (cnt%viz_freq == 0):
            ss.plot.figure(1)
            plt.clf()
            ss.plot.plot(pvar,'b.-')
            ss.plot.plot('mlbeta','k--')
            
            ss.plot.figure(2)
            plt.clf()
            ss.plot.plot('u','b.-')

ML_data_trim = ML_data[:,:mlcnt]
            
if teach:
    numpy.savetxt("bulk_sod.txt",ML_data_trim.T,fmt="%.4e",newline="\n",delimiter=" ")
        
ss.writeGrid()
ss.write()
