# Run a bunch of test and check answers
import os,sys
import numpy as npy
import subprocess
from testObj import *


test_dir = os.path.dirname(os.path.abspath(__file__))
bin_dir  = os.path.join( test_dir, '../bin')
root_dir = os.path.join( test_dir, '..')
pyranda_exe = sys.executable

tests = []    # List of test objects
dbase = {}    # Dictionary of baselines
relE  = {}    # Dictionary of relative errors for baselines

# Add tests here
execfile('cases/test1DAdvection.py')
execfile('cases/testIntegration.py')
execfile('cases/testMM_simple.py')
execfile('cases/test2deuler.py')
execfile('cases/testCylinder.py')
execfile('cases/testHeat1D.py')
execfile('cases/testTaylorGreen.py')


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
        baseline = dbase[test.name]  # Baseline value/file
        try:
            rdiff = relE[test.name]  # Rel. error for this baseline
        except:
            rdiff = 1.0e-4
            
        
        if curve:
            # Check curve
            diff = checkProfile( baseline, pout)
        else:
            # Check if scalar compare    
            diff = checkScalar( float(baseline) , pout)

        if diff < rdiff:
            testPass = True
            print 'Pass: (Rel. Error = %s )' % diff
            fout = '%s -- %s' % (test.name,pout)
            print fout
            passed += 1
            if curve:
                sexe("rm %s" % pout,ret_output=False,echo=False)
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

if failed > 0:
    print '\n\n\n===== New baselines ====='
    print new_baselines

plt.show()
    
        
