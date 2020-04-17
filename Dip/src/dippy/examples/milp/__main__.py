'''
Created on Dec 29, 2013

@author: ted
'''
from __future__ import division
from builtins import str
from past.utils import old_div
from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize
try:
    from src.dippy import Solve
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import Solve
    from coinor.dippy.examples.gen_func import *

from .milp_func import *

args = parseArgs()

tol = pow(pow(2, -24), old_div(2.0, 3.0))

prob = formulate(args)

dippyOpts = addDippyOpts(args)

dippyOpts['TolZero'] = str(tol)
#dippyOpts['LogDebugLevel'] = '3'
#dippyOpts['LogLevel'] = '4'
#dippyOpts['LogDumpModel'] = '5'
#dippyOpts['ALPS'] = {'msgLevel' : 3}
    
Solve(prob, dippyOpts)

if prob.display_mode != 'off':
    prob.Tree.display()

