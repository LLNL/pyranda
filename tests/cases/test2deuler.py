# Copy and paste baselines here
baselines = """
euler-2d-64 -- euler-2d-64.dat
euler-2d-128 -- euler-2d-128.dat
KH-2d-64 -- KelvinHelmholtzKH-2d-64.dat
"""

# Update dictionary of baseline scalars
dbase.update( baseDict( baselines) )

test = testObj('euler-2d-64')
test.script = 'examples/euler.py'
test.args = ['64','1',test.name]
#tests.append( test )

test = testObj('euler-2d-128')
test.script = 'examples/euler.py'
test.args = ['128','1',test.name]
#tests.append( test )

# Kelvin-Helmholtz test problem
test = testObj('KH-2d-64')
test.script = 'examples/KH.py'
test.args = ['64','0',test.name,"1"]
tests.append( test )

