#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from past.utils import old_div
from .bin_pack_func import BinPackProb, formulate, solve

import sys

if __name__ == '__main__':
    # Python starts here
    bpp = BinPackProb(ITEMS  = [1, 2, 3, 4, 5, 6],
                      volume = {1: 2, 2: 5, 3: 3, 4: 3, 5: 3, 6: 2},
                      capacity = 9)
  
    prob = formulate(bpp)
  
    prob.tol = pow(pow(2, -24), old_div(2.0, 3.0))
    if len(sys.argv) > 1:
        xopt = solve(prob, sys.argv[1])
    else:
        xopt = solve(prob)
  
    if xopt is not None:
        for var in prob.variables():
            print(var.name, "=", xopt[var])
    else:
        print("Dippy could not find and optimal solution")
  
    if prob.display_mode != 'off':
        numNodes = len(prob.Tree.get_node_list())
        if (prob.Tree.attr['display'] == 'pygame') or (prob.Tree.attr['display'] == 'xdot'):
            prob.Tree.display()
