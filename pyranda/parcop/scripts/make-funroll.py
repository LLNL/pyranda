#!/bin/local/bin/python
import os,sys
import shutil
import glob


# Select files for code-to-code transforms
pwd = os.path.dirname(os.path.abspath(__file__))
DIRS = []
DIRS.append(pwd + '/../../parcop')  ### Add source file dirs here ###


#### Addd your Fortran suffix support here
fsuff = ".f90"

try:
    option = int(sys.argv[1])
except:
    print("Give a valid option (1 or 2)")

try:
    fileWildCard = sys.argv[2]
except:
    fileWildCard = None


# Option-1: Make conversion (safe to call multple times)
Files = []
for D in DIRS:
    Files += glob.glob( D + '/*%s' % fsuff )
if fileWildCard:
    Files = [ F for F in Files if fileWildCard in F]

if option == 1:
    # Not go through each and convert it
    for ff in Files:
        try:
            os.system('python %s/funroller.py %s %s %s' % (pwd,
                                                        ff,
                                                        fsuff,
                                                        0) )
        except:
            print("Error processing file:%s" % ff)
    exit()


# Option-2: Revert space (safe to call multple times)
Files = []
for D in DIRS:
    Files += glob.glob( D + '/*.fexl' )
if fileWildCard:
    Files = [ F for F in Files if fileWildCard in F]
            
# Option-2: Revert space (safe to call multple times)
if option == 2:    
    for ff in Files:
        cmd = 'mv -f %s %s' % (ff,ff.replace(fsuff+".fexl",fsuff))
	print("Reverting %s" % ff)
        os.system( cmd )
    exit()


        
