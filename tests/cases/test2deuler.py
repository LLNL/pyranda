# Copy and paste baselines here
baselines = """
euler-2d-64 -- euler-2d-64.dat
euler-2d-128 -- euler-2d-128.dat
"""

# Update dictionary of baseline scalars
dbase.update( baseDict( baselines) )

test = testObj('euler-2d-64')
test.script = 'examples/euler.py'
test.args = ['64','1',test.name]
tests.append( test )

test = testObj('euler-2d-128')
test.script = 'examples/euler.py'
test.args = ['128','1',test.name]
tests.append( test )




