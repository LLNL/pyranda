# Run a bunch of test and check answers
import os,sys
import numpy as npy
import subprocess
import matplotlib.pyplot as plt

def sexe(cmd,ret_output=False,echo = False):
    """ Helper for executing shell commands. """
    if echo:
        print "[exe: %s]" % cmd
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        res =p.communicate()[0]
        return p.returncode,res
    else:
        return subprocess.call(cmd,shell=True)


def baseDict(string):
    dbase = {}
    # Make dictionary
    for bb in string.split('\n'):
        if bb:
            name = bb.split('--')[0].strip()
            diff = bb.split('--')[1].strip()
            dbase[name] = diff
    return dbase
            
def checkScalar(baseline,pout):
    if abs(baseline) > 0.0:
        diff = npy.abs( float(pout)-baseline ) / npy.abs( baseline )
    else:
        diff = npy.abs( float(pout)-baseline )
    return diff


def checkProfile( baseline, pout ):

    baseline = 'baselines/%s' % baseline

    base = npy.loadtxt( baseline )
    test = npy.loadtxt( pout )

    diff = npy.abs( base[1,:] - test[1,:] ) / npy.max( npy.abs( base[1,:] ) )

    m_diff = npy.max( diff )

        
    return m_diff

def plotError( baseline, pout ):

    baseline = 'baselines/%s' % baseline

    base = npy.loadtxt( baseline )
    test = npy.loadtxt( pout )

    diff = npy.abs( base[1,:] - test[1,:] ) / npy.max( npy.abs( base[1,:] ) )
    
    plt.figure()
    plt.plot([1,2,3])
    plt.subplot(121)
    plt.plot(base[0,:],base[1,:],'k-')
    plt.plot(test[0,:],test[1,:],'b-')

    plt.subplot(122)
    plt.plot(base[0,:],diff,'r-')

    
class testObj:

    def __init__(self,name):

        self.name = name
        self.parallel = False
        self.np = 1
        self.script = ''
        self.args = None
