# Copy and paste baselines here
baselines = """
shu-osher-CR -- 1.57446677695
"""

# Update dictionary of baseline scalars
dbase.update( baseDict( baselines) )

test = testObj('shu-osher-CR')
test.script = 'examples/shu_osher.py'
test.args = ['1']
tests.append( test )




