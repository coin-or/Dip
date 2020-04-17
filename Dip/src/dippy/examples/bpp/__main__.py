#!/usr/bin/env python

# bin_pack_instance.py
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from past.utils import old_div
from .bin_pack_func import BinPackProb, formulate, solve
from os.path import dirname
from inspect import getfile
import coinor.dippy.examples.bpp
import importlib as ilib

import sys

if len(sys.argv) > 1:
    if sys.argv[1] == '-h' or sys.argv[1] == '--help' or len(sys.argv) > 2:
        print('Usage: bpp <module_name>')
        print('       module_name : Python module containing instance data')
        print('                     For example file, check directory')
        print('                    ', dirname(getfile(coinor.dippy.examples.bpp)))
        exit()
    else:
        m = ilib.import_module(sys.argv[1])
        bpp = BinPackProb(ITEMS = m.ITEMS, volume = m.volume, capacity = m.capacity)
        prob = formulate(bpp)
else:
    bpp = BinPackProb(ITEMS  = [1, 2, 3, 4, 5],
                      volume = {1: 2, 2: 5, 3: 3, 4: 7, 5: 2},
                      capacity = 8)
    prob = formulate(bpp)

xopt = solve(prob)
  
if xopt is not None:
    for var in prob.variables():
        print(var.name, "=", xopt[var])
else:
    print("Dippy could not find and optimal solution")    
  
if prob.display_mode != 'off':
    numNodes = len(prob.Tree.get_node_list())
    if ((prob.Tree.attr['display'] == 'pygame') or 
        (prob.Tree.attr['display'] == 'xdot')):
        prob.Tree.display()


