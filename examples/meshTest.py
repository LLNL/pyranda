import numpy
import matplotlib.pyplot as plt




def zoomMesh(npts,x1,xn,xa,xb,tt,dxc,dxf):

    x = numpy.zeros( (npts) )

    x[0] = x1
    dx = dxc
    for i in range(1,npts):
        
        x[i] = x[i-1] + dx

        xmin = (xa+xb)/2.0
        xh = (xb-xa)/2.0
        xpr  = numpy.abs(x[i]-xmin)
        w = 0.5*(numpy.tanh( (xpr - xh)/tt ) + 1.0)
        dx = w*dxc + (1.0-w)*dxf

    return x


def expMesh(npts,x1,xn,xa,xb,dxf,iexp):

    BmL = ( numpy.log( (xa-x1)/dxf  * iexp + 1.) / iexp )
    Bint = numpy.floor( BmL )

    nC = (xb-xa)/dxf
    
    x = numpy.zeros( (npts) )
    x[0] = x1

    xAA = xa
    for i in range(1,npts):

        if i <= Bint:
            dx = dxf*numpy.exp( (Bint-(i-1))*iexp )
            xAA = x[i-1] + dx 
        elif (i>Bint) and (i<Bint+nC):
            dx = dxf
        elif (i>=Bint+nC):
            dx = dxf*numpy.exp( (i-(Bint+nC))*iexp )
                    
        x[i] = x[i-1] + dx        
        
    x += (xa-xAA)    
    return x
        
    
def efunc(npts,x1,xn,xa,xb,dxf,iexp):
    xx = expMesh(npts,x1,xn,xa,xb,dxf,iexp)
    error = xx[-1] - xn
    return error


def sfunc(npts,x1,xn,xa,xb,tt,dxc,dxf):
    xx = zoomMesh(npts,x1,xn,xa,xb,tt,dxc,dxf)
    error = xx[-1] - xn
    return error
    

def zoomMesh_solve(npts,x1,xn,xa,xb,tt,dxf):

    # take initial guess at dxc and solve
    nmin = (xb-xa)/dxf
    if nmin > npts:
        print("Not enough points given")
        exit()
    dxc = ((xn-x1) - (xb-xa)) / (npts-nmin)

    # NR solvers
    delta = 1.0001
    f1 = xn
    cnt = 1
    while ( (abs(f1) > .01*xn) and cnt<100):
        f1 = sfunc(npts,x1,xn,xa,xb,tt,dxc,dxf)
        f2 = sfunc(npts,x1,xn,xa,xb,tt,dxc*delta,dxf)
        #f1 = efunc(npts,x1,xn,xa,xb,dxf,exp)
        #f2 = efunc(npts,x1,xn,xa,xb,dxf,exp*delta)        
        dxc = dxc*(1.0- f1/(f2-f1)*(delta-1.0))
        cnt += 1
    
        #print(zoomMesh(npts,x1,xn,xa,xb,tt,dxc,dxf)[-1])
    return zoomMesh(npts,x1,xn,xa,xb,tt,dxc,dxf)


def expMesh_solve(npts,x1,xn,xa,xb,dxf,exp):

    # take initial guess at dxc and solve
    nmin = (xb-xa)/dxf
    if nmin > npts:
        print("Not enough points given")
        exit()
    dxc = ((xn-x1) - (xb-xa)) / (npts-nmin)

    # NR solvers
    delta = 1.0001
    f1 = xn
    cnt = 1
    while ( (abs(f1) > .01*xn) and cnt<100):
        f1 = efunc(npts,x1,xn,xa,xb,dxf,exp)
        f2 = efunc(npts,x1,xn,xa,xb,dxf,exp*delta)        
        exp = exp*(1.0- f1/(f2-f1)*(delta-1.0))
        cnt += 1
    
        print(f1,f2)
        #print(zoomMesh(npts,x1,xn,xa,xb,tt,dxc,dxf)[-1])
    return expMesh(npts,x1,xn,xa,xb,dxf,exp)




npts = 100
x1 = 0
xn = 1
xa = .4
xb = .6
tt = .02
dxf = (xn-x1) / float(npts) * .75

xx = zoomMesh_solve(npts,x1,xn,xa,xb,tt,dxf)
plt.plot( xx[:-1], numpy.diff(xx) )
plt.figure()
plt.plot( xx )
plt.show()

npts = 100
x1 = -1
xn = 1
xa = .4
xb = .6
dxf = (xb-xa) / float(npts) * 1.3
exp = 0.3

xx = expMesh_solve(npts,x1,xn,xa,xb,dxf,exp)
#xx = expMesh(npts,x1,xn,xa,xb,dxf,exp) 
plt.plot( xx[:-1], numpy.diff(xx) ,'b-o')
plt.figure()
plt.plot( xx ,'b-o')
plt.show()



