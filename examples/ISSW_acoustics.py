from __future__ import print_function
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm


from pyranda import pyrandaSim, pyrandaBC, pyrandaTimestep, pyrandaIBM


# Try to get args for testing
try:
    Npts = int(sys.argv[1])
except:
    Npts = 64

try:
    test = bool(int(sys.argv[2]))
except:
    test = False

try:
    testName = (sys.argv[3])
except:
    testName = None


## Define a mesh
#Npts = 32
L = numpy.pi * 2.0  
dim = 1
gamma = 1.4
problem = '1D_LinearReflection'



rho0 = 1.0
p0   = 1.0
gamma = 1.4
mach = 2.0
s0 = numpy.sqrt( p0 / rho0 * gamma )
u0 = s0 * mach
e0 = p0/(gamma-1.0) + rho0*.5*u0*u0
k0 = 0.1     # Wave number of sinewave
amp = 1.0e-4

# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -div(:rho:*:u:) 
ddt(:rhou:) =  -div(:rhou:*:u: + :p: - :tau:) 
ddt(:Et:)   =  -div( (:Et: + :p: - :tau:)*:u: - :tx:*:kappa: )
# Level set equation
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity (old school way)
:div:       =  div(:u:)
:beta:      =  gbar( ring(:div:) * :rho: ) * 7.0e-3
:tau:       =  :beta: * :div: 
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
# Apply constant BCs
[:u:,:v:,:w:] = ibmV( [:u:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:], [:u1:,:u2:,0.0] )
:rho: = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
:p:   = ibmS( :p:   , :phi:, [:gx:,:gy:,:gz:] )
:right: = simtime-( ((2.0*k0+0.5)*2.0*pi + Xshift + meshx)/s0)
:left:  = simtime-( ( 2.0*xWall + (2.0*k0+0.5)*2.0*pi + Xshift - meshx) / s0)
:inlet: = p0*amp*exp( -:right:**2 / (k0*pi/s0)**2 ) 
:inlet: = :inlet: + p0*amp*exp( -:left:**2 / (k0*pi/s0)**2 ) 
:pA: = :inlet: + p0
:error: = abs(:p: - :pA:)
bc.extrap(['rho','p','u'],['xn'])
bc.field(['p'],  ['x1'],p0 + :inlet:)
bc.field(['rho'],['x1'],rho0 + :inlet:/(s0*s0))
bc.field(['u'],  ['x1'],:inlet:/(s0*rho0))
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: )
:rhou: = :rho:*:u:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.2* dt.diff(:beta:,:rho:)
:dt: = numpy.minimum(:dt:,:dtB:)
:umag: = sqrt( :u:*:u: )
"""



# Initialize variables
ic = """
:gamma: = 1.4
:R: = 1.0
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
:phi: = xWall - meshx
:rho: = 1.0 + 3d()
:p:  =  1.0 + 3d() 
:u: = 0.0 + 3d() 
:Et: = :p:/( :gamma: - 1.0 ) + .5*:rho:*(:u:*:u: )
:rhou: = :rho:*:u:
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
[:gx:,:gy:,:gz:] = grad( :phi: )
"""

def solve(Npts,test=True,delta=0.0):


    Lp = L * (Npts-1.0) / Npts

    mesh_options = {}
    mesh_options['coordsys'] = 3
    mesh_options['periodic'] = numpy.array([False, False, True])
    mesh_options['dim'] = 3
    mesh_options['x1'] = [ 0.0 , 0.0  ,  0.0 ]
    mesh_options['xn'] = [ Lp  ,  Lp  ,  Lp ]
    mesh_options['nn'] = [ Npts, 1    ,  1  ]

    delta = delta * (Lp/Npts)

    dx = (L/Npts)
    try:
        myOff = Offset
    except:
        myOff = .5
    xWall = (int(.75*Npts)+myOff)*dx
    
    Xshift = -2*(xWall - .75*2*numpy.pi)
    
    # Initialize a simulation object on a mesh
    ss = pyrandaSim(problem,mesh_options,silent=test)
    ss.addPackage( pyrandaBC(ss) )
    ss.addPackage( pyrandaIBM(ss) )
    ss.addPackage( pyrandaTimestep(ss) )

    # Add the EOM to the solver
    meom = eom.replace('u0',str(u0)).replace('p0',str(p0)).replace('rho0',str(rho0)).replace('s0',str(s0)).replace('k0',str(k0)).replace('xWall',str(xWall)).replace('amp',str(amp))
    meom = meom.replace('Xshift',str(Xshift))
    ss.EOM(meom)

    mic = ic.replace('mach',str(mach)).replace('xWall',str(xWall))

    # Set the initial conditions
    ss.setIC(mic)

    # Write a time loop
    time = 0.0
    viz = True

    # Approx a max dt and stopping time
    tt = 10.0 

    # Start time loop
    cnt = 1
    viz_freq = 1#int(Npts/5)
    pvar = 'p'


    CFL = 1.0
    dt = ss.variables['dt'].data * CFL*.01

    #ss.parse(":p: = sin(meshx*2.0)")
    #ss.parse(":p:   = ibmS2( :p:   , :phi:, [:gx:,:gy:,:gz:] )")
    #ss.plot.plot('p')
    
    #import pdb
    #pdb.set_trace()
    
    error = 0.0
    dtSample = 1.0e5
    while tt > time:

        # Update the EOM and get next dt
        time = ss.rk4(time,dt)
        dt = min( ss.variables['dt'].data * CFL, dt*1.1)
        dt = min(dt, (tt - time) )


        #if time > 9.0 and time < 11.0:
        #error = max( error, numpy.abs( 1.01 - ss.eval("max(:p:)") ) )
            

        
        # Print some output
        ss.iprint("%s -- %s --- %f" % (cnt,time,dt)  )
        cnt += 1
        if viz and (not test):
            v = ss.PyMPI.zbar( ss.variables[pvar].data )
            phi = ss.PyMPI.zbar( ss.variables['phi'].data )

            if (cnt%viz_freq == 0) :#or True:
                ss.plot.figure(2)
                ss.plot.clf()
                plt.ylim( (1.0,1.022) )
                ss.plot.plot('p','ko')
                ss.plot.plot('pA','b--')
                xs = xWall #3.*numpy.pi/2.
                plt.plot( [xs,xs],[0,1.02],'k--')


                ss.plot.figure(1)
                ss.plot.clf()
                ss.plot.plot('u')
                plt.plot( [xs,xs],[-.01,.01],'k--')


                #ss.plot.figure(3)
                #ss.plot.clf()
                #ss.plot.plot('rho')

                #ss.plot.figure(4)
                #ss.plot.clf()
                #ss.plot.plot('Et')


                
                plt.pause(.01)
                try:
                    input("HERE")
                except:
                    pass
                
    error = ss.eval("mean( where(:phi:>0,:error:,0.0))")/amp
    minPhi = ss.eval("min(where(:phi:>0,:phi:,1e6)    )")
    print(error,minPhi)

    #if (Npts == 800):
    #    ss.plot.plot('p','k-',linewidth=2.5,label='$N=800$')
    #if (Npts == 400):
    #    ss.plot.plot('p','k--',linewidth=2.5,label='$N=400$')
    
    return [ss,error,minPhi]
    

RES = [25,50,100,200,400,800]

#RES = [50,51,52,53,54,55,56,57,58,59,60]

errors = []
mySS = None
for r in RES:
    [ss,ee,mphi] = solve(r,test=True)
    errors.append( ee )
    if (r == 200):
        mySS = ss

errors = numpy.array(errors)
RES = numpy.array(RES)
CR = numpy.log( (errors[1:]/errors[:-1])) / numpy.log( (RES[1:]/RES[:-1]))
print(CR)

DX = 1.0 / numpy.array( RES) * 2.0*numpy.pi


#DEL = [.0,.1,.2,.3,.4,.5,.6,.7,.8,.9]

#for delta in DEL:
#    [ss,ee,mphi] = solve(60,test=True,delta=delta)

#[ss,ee,mphi] = solve(50,test=False)

#runMe = """

LW = 2.5
from matplotlib import rcParams
import os
rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16
rcParams['legend.fontsize'] = 20
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
floc = "/Users/olson45/Documents/Conference Travel/ISSW32_2019/ISSW32_bjo/figures"


#Offset = .1
#execfile("ISSW_acoustics.py")
#plt.loglog(1./DX,errors,'bo-', linewidth=LW)

Offset = .5
#execfile("ISSW_acoustics.py")
plt.figure(2)
plt.loglog(1./DX,errors,'ko-', linewidth=LW)

#Offset = .9
#execfile("ISSW_acoustics.py")
#plt.loglog(1./DX,errors,'ro-', linewidth=LW)


rate = 2
scl = 1./DX[0]**(-rate)/errors[0]
plt.loglog( 1./DX , .5/scl * (1./DX)**(-rate), 'k--')

plt.xlabel(r"$\frac{1}{\Delta x}$") 
plt.ylabel(r"$L_2$ Error")
plt.tight_layout(pad=1.0)
plt.savefig(os.path.join(floc,'acoust_conv.pdf'))


plt.figure(1)

mySS.plot.plot('p','bo',linewidth=2.5)
mySS.plot.plot('pA','k--',linewidth=2.5)
plt.xlim((1.0,3.0))

plt.xlabel(r"$x$") 
plt.ylabel(r"$p(x)$")
plt.tight_layout(pad=1.0)
plt.savefig(os.path.join(floc,'acoust_prof.pdf'))


plt.show()


#"""
