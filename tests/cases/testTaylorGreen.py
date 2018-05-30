# Copy and paste baselines here
baselines = """
taylorgreen -- 1.00127724838
"""

# Update dictionary of baseline scalars
dbase.update( baseDict( baselines) )

script = 'examples/TaylorGreen.py'

test = testObj('taylorgreen')
test.script = script
test.args = ['32','1']
tests.append( test )

