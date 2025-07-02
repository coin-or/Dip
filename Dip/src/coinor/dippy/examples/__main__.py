from os.path import dirname
from inspect import getfile
import coinor.dippy.examples

print('''
Welcome the DiPPy examples! Here you will find a wide range of examples.
The examples are installed in the directory
''')
print(dirname(getfile(coinor.dippy.examples)))
print(
'''
They're meant to be inspected, so please do look at the source code!
It's also available in the source repo

https://github.com/coin-or/Dip/tree/master/Dip/src/dippy/examples

Each example can be run directly as a Python module with a default input data, e.g.,

python -m coinor.dippy.examples.cflp

For help on how to run a given example, run the module with '--help'.
Most examples also accept input data in the form of a Python module, e.g.,

python -m coinor.dippy.examples.cflp coinor.dippy.examples.bpp.bpp_data

The example modules available are

bpp:     Bin Packing Problem
cflp:    Capacitated Facility Location Problem
coke:    Steel Manufacturing
csp:     Cutting Stock Problem
cvpmp:   Capacitated P-median
gap:     Generalized Assignment
milp:    Random Block Structured MILP
tsp:     Traveling Salesman Problem
wedding: Wedding Planner

Enjoy!
''')
