'''
Created on Dec 29, 2013

@author: ted
'''
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div

#!/usr/bin/env python

import sys
from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize

try:
    from src.dippy import Solve
except ImportError:
    from coinor.dippy import Solve

from .milp_func import *
    
algo = 'Cut'
rand_seed = 2
if len(sys.argv) >1:
    if sys.argv[1] == '-h' or sys.argv[1] == '--help' or len(sys.argv) < 5 or len(sys.argv) > 7:
        print('Usage: milp [ i ] [ j ] [ k ] [ l] [ algo ] [seed]')
        print('       i : number of blocks')
        print('       j : number of variables per block')
        print('       k : number of constraints per block')
        print('       l : number of linking constraints')
        print("       algo: Algorithm 'Cut', 'Price', 'PriceCut'")
        print('       seed: random seed')
        exit()
    elif len(sys.argv) > 1:
        numBlocks = int(sys.argv[1])
        numVarsPerBlock = int(sys.argv[2])
        numConsPerBlock = int(sys.argv[3])
        numLinkingCons = int(sys.argv[4])
        if len(sys.argv) >= 6:
            algo = sys.argv[5]
        if len(sys.argv) == 7:
            rand_seed = sys.argv[6]
else:
    numBlocks = 10
    numVarsPerBlock = 10
    numConsPerBlock = 5
    numLinkingCons = 10

tol = pow(pow(2, -24), old_div(2.0, 3.0))

prob = formulate(numBlocks, numVarsPerBlock, numConsPerBlock, numLinkingCons, rand_seed)

dippyOpts = {}
if algo == 'PriceCut':
    dippyOpts['doPriceCut'] = '1'
    dippyOpts['CutCGL'] = '1'
elif algo == 'Price':
    dippyOpts['doPriceCut'] = '1'
    dippyOpts['CutCGL'] = '0'
else:
    dippyOpts['doCut'] = '1'

dippyOpts['TolZero'] = str(tol)
dippyOpts['SolveMasterAsIp'] = '0'
dippyOpts['generateInitVars'] = '1'
#dippyOpts['LogDebugLevel'] = '3'
#dippyOpts['LogLevel'] = '4'
#dippyOpts['LogDumpModel'] = '5'
#dippyOpts['ALPS'] = {'nodeLogInterval' : 1,
#                     'nodeLimit' : 1,
#                     'msgLevel' : 3
#                     }
    
Solve(prob, dippyOpts)

if prob.display_mode != 'off':
    prob.Tree.display()

