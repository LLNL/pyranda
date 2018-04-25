# Run a bunch of test and check answers
import os,sys
import numpy as npy
import subprocess
from testObj import *


test_dir = os.path.dirname(os.path.abspath(__file__))
bin_dir  = os.path.join( test_dir, '../bin')
root_dir = os.path.join( test_dir, '..')
pyranda_exe = 'python' #os.path.join( bin_dir, 'pyranda')
pyranda_exe = '/opt/local/bin/python2.7'
pyranda_mpi = os.path.join( bin_dir, 'pympirun')



tests = []    # List of test objects
dbase = {}    # Dictionary of baselines

# Add tests here
execfile('cases/test1DAdvection.py')
execfile('cases/testMM_simple.py')
execfile('cases/test2deuler.py')
execfile('cases/testCylinder.py')
execfile('cases/testHeat1D.py')


summary = ''
passed = 0
failed = 0
new_baselines = ''


# Run tests
for test in tests:

    script = os.path.join(root_dir,test.script)

    exe = pyranda_exe
    if test.parallel:
        exe = pyranda_mpi + ' -n %s %s' % (test.np,exe)
    
    # Args
    sargs = ''
    for arg in test.args:
        sargs += '%s ' % arg

        
    cmd = '%s %s %s' % (exe,script,sargs)

    out = sexe(cmd,ret_output=True,echo=False)
    pout = out[1].split('\n')[-2]
    curve = False
    if '.dat' in pout:
        curve = True
    
    # Diff against baseline
    try:
        baseline = dbase[test.name]

        if curve:
            # Check curve
            diff = checkProfile( baseline, pout)
        else:
            # Check if scalar compare    
            diff = checkScalar( float(baseline) , pout)

        if diff < 1.0e-4:
            testPass = True
            print 'Pass: (Rel. Error = %s )' % diff
            fout = '%s -- %s' % (test.name,pout)
            print fout
            passed += 1
        else:
            testPass = False
            print 'Fail: (Rel. Error = %s )' % diff
            fout = '%s -- %s' % (test.name,pout)
            print fout
            new_baselines += fout + '\n'
            failed += 1

            if curve:
                plotError( baseline, pout )
            
            

    except:
        testPass = False
        print 'Fail: (No baseline data found )'
        fout = '%s -- %s' % (test.name,pout)
        print fout
        new_baselines += fout + '\n'
        failed += 1

print "\n\n\n=====Testing Summary======"
print "Passed: %s" % passed
print "Failed: %s" % failed

print '\n\n\n===== New baselines ====='
print new_baselines

plt.show()

    
        
