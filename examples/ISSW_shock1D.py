import re
import sys
import time
import numpy 
import matplotlib.pyplot as plt
from matplotlib import cm

from pyranda import pyrandaSim, pyrandaBC,pyrandaIBM,pyrandaTimestep
from pyranda.pyranda import pyrandaRestart


## Define a mesh
L = numpy.pi * 2.0
Npts = 100


Ms = 1.2
gamma = 1.4
p0 = 1.0
rho0 = 1.0
c0 = numpy.sqrt( p0 * gamma  / rho0 )
Et0 = p0/(gamma-1.0) 

gm1 = (gamma-1.0)
gp1 = (gamma+1.0)
pRatio = (2.*gamma*Ms**2-gm1) / gp1
rhoRatio = gp1*Ms**2/( gm1*Ms**2+2.0)
uRatio = (gm1*Ms**2+2.)/(gp1*Ms**2)

p1 = p0 * pRatio
rho1 = rho0 * rhoRatio
c1 = numpy.sqrt( p1 * gamma / rho1 )
M1 = numpy.sqrt( (gm1*Ms**2 + 2.0)/( 2.*gamma*Ms**2 - gm1 ) )

u1 = numpy.abs( Ms*c0 - M1*c1 )
Et1 = p1/(gamma-1.0) + .5*rho1*u1*u1


# Define the equations of motion
eom ="""
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: - :tx:*:kappa:)
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
:T:         = :p: / (:rho: * :R: )
# Artificial bulk viscosity
:div:       =  ddx(:u:) 
:beta:      =  gbar( ring(:div:) * :rho:) * 7.0e-2
:tau:       =  :beta:*:div:
[:tx:,:ty:,:tz:] = grad(:T:)
:kappa:     = gbar( ring(:T:)* :rho:*:cv:/(:T: * :dt: ) ) * 1.0e-3
# Apply the IBM
[:u:,:v:,:w:] = ibmV( [:u:,:v:,0.0], :phi:, [:gx:,:gy:,:gz:])
:rho: = ibmS( :rho: , :phi:, [:gx:,:gy:,:gz:] )
:p:   = ibmS( :p:   , :phi:, [:gx:,:gy:,:gz:] )
:Et:  = :p: / ( :gamma: - 1.0 )  + .5*:rho:*(:u:*:u: + :v:*:v:)
:rhou: = :rho:*:u:
:rhov: = :rho:*:v:
# Apply constant BCs
bc.const(['Et']  ,['x1'],Et1    )
bc.const(['rho'] ,['x1'],rho1   )
bc.const(['rhou'],['x1'],rho1*u1)
bc.const(['p']   ,['x1'],p1     )
# Time step control
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)
:dtB: = 0.2* dt.diff(:beta:,:rho:)
:dt: = numpy.minimum(:dt:,:dtB:)
"""



# Initial conditions 
ic = """
:R: = 1.0
:cp: = :R: / (1.0 - 1.0/:gamma: )
:cv: = :cp: - :R:
:phi: = xWall - meshx
[:gx:,:gy:,:gz:] = grad( :phi: )
:gamma: = myGamma
:rho:  = gbar( where( meshx < pi, rho1    , rho0  ))
:rhou: = gbar( where( meshx < pi, rho1*u1 , 0.0   ))
:Et:   = gbar( where( meshx < pi, Et1     , Et0   ))
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*.1
"""


xWall = .75*L 

def solve(Npts,test=True,delta=0.5):


    Lp = L * (Npts-1.0) / Npts

    imesh = """
    xdom = (0.0, Lp, Npts)
    """.replace('Lp',str(Lp)).replace('Npts',str(Npts))


    # Initialize a simulation object on a mesh
    ss = pyrandaSim('1D_shock_reflection',imesh)
    ss.addPackage( pyrandaBC(ss) )
    ss.addPackage( pyrandaIBM(ss) )
    ss.addPackage( pyrandaTimestep(ss) )


    dx = (L/Npts)
    try:
        myOff = delta
    except:
        myOff = .5

    xWall = .75*L #(int(.75*Npts)+myOff)*dx


    # Add the EOM to the solver
    myEOM = eom.replace('u1',str(u1)).replace('rho1',str(rho1))
    myEOM = myEOM.replace('Et1',str(Et1)).replace('p1',str(p1))
    ss.EOM(myEOM)

    # Set the initial conditions
    myIC =   ic.replace('rho1',str(rho1)).replace('rho0',str(rho0))
    myIC = myIC.replace('Et1',str(Et1)).replace('Et0',str(Et0))
    myIC = myIC.replace('u1',str(u1)).replace('myGamma',str(gamma))
    myIC = myIC.replace('xWall',str(xWall))
    ss.setIC(myIC)



    # Write a time loop
    time = 0.0
    tt = 2.0

    # Start time loop
    dt = ss.var('dt').data
    cnt = 1
    viz_freq = 10
    pvar = 'rho'
    viz = True

    while tt > time:

        # Update the EOM and get next dt
        time = ss.rk4(time,dt)
        dt = min( .5*ss.var('dt').data, (tt - time) )

        # Print some output
        ss.iprint("%s -- %s" % (cnt,time)  )
        cnt += 1
        if viz and (not test):

            if (cnt%viz_freq == 0):
                ss.plot.figure(1)
                plt.clf()
                xs = xWall #3.*numpy.pi/2.
                plt.plot( [xs,xs],[-0.1,2.35],'k--')
                ss.plot.plot('rho','b.-')
                ss.plot.plot('p','r.-')
                #ss.plot.figure(2)
                #plt.clf()
                ss.plot.plot('u','k--')


    #ss.writeGrid()
    #ss.write()
    return ss


RES = [50,100,200,400,800,1600]
sols = []
for r in RES:
    ss = solve(r,test=True)
    sols.append( ss )


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

    
# Plot all the solutions
xs = xWall
for ss in sols:
    ss.plot.plot('p','k.-',linewidth=LW)
    plt.plot( [xs,xs],[1.45,2.35],'k--')

plt.xlim( (3.0,5.3) )
plt.ylim( (1.45,2.35) )
plt.xlabel(r"$x$") 
plt.ylabel(r"$p(x)$")
plt.tight_layout(pad=1.0)
plt.savefig(os.path.join(floc,'shock1D_p.pdf'))


#plt.show()
# Compute error from ref.
ref = sols[-1]
uref = ref.var('u').data[:,0,0]
pref = ref.var('p').data[:,0,0]
xref = ref.x().data[:,0,0]
pL2s = []
uL2s = []
DX = []
for ss in sols[:-1]:
    p = ss.var('p').data[:,0,0]
    u = ss.var('u').data[:,0,0]
    x = ss.x().data[:,0,0]
    pp = numpy.interp(x,xref,pref)
    uu = numpy.interp(x,xref,uref)
    mask = numpy.where( x>3.0, 1.0, 0.0)
    mask = numpy.where( x<xWall, mask, 0.0)
    errorP = numpy.abs( p - pp ) * mask
    errorU = numpy.abs( u - uu ) * mask
    pL2 = numpy.mean( errorP )
    pLi = numpy.max( errorP )
    uL2 = numpy.mean( errorU )
    uLi = numpy.max( errorU )
    print("L2: %s ---- Li: %s" % (pL2,pLi) )
    print("L2: %s ---- Li: %s" % (uL2,uLi) )
    pL2s.append( pL2 )
    uL2s.append( uL2 )
    DX.append(  ss.dx )
    

DX   = numpy.array( DX   )
pL2s = numpy.array( pL2s )

Lref = xWall - 3.7

plt.figure(2)
plt.loglog( Lref / DX, pL2s,'ko-', linewidth=LW)

rate = 1
scl = 1./DX[0]**(-rate)/pL2s[0]
plt.loglog( 1./DX , .5/scl * (1./DX)**(-rate), 'k--')

plt.xlabel(r"$\frac{L_w}{\Delta x}$") 
plt.ylabel(r"$L_2$ Error")
plt.tight_layout(pad=1.0)
plt.savefig(os.path.join(floc,'shock1D_conv.pdf'))
