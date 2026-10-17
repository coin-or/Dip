from builtins import range
from builtins import object
import sys

from pulp import LpVariable, lpSum, LpBinary, LpStatusOptimal

try:
    import path
except ImportError:
    pass
        
try:
    import dippy
except ImportError:
    try:
        import src.dippy as dippy
    except ImportError:
        import coinor.dippy as dippy

from math import floor, ceil

class MDBinPackProb(object):
    def __init__(self, ITEMS, LIMITS, volume, capacity):
        self.ITEMS = ITEMS
        self.LIMITS = LIMITS
        self.volume = volume
        self.capacity = capacity

        self.BINS = list(range(len(ITEMS))) # Create 1 bin for each item, indices 
                                      # start at 0
    
def formulate(bpp):
    prob = dippy.DipProblem("Bin Packing",
                            display_mode = 'xdot',
#                           layout = 'bak',
                            display_interval = None,
                            )

    assign_vars = LpVariable.dicts("x",
                                   [(i, j) for i in bpp.ITEMS
                                    for j in bpp.BINS],
                                   cat=LpBinary)
    use_vars    = LpVariable.dicts("y", bpp.BINS, cat=LpBinary)
    waste_vars  = LpVariable.dicts("w", [(j, k) for j in bpp.BINS
                                          for k in bpp.LIMITS], 0, None)

    prob += lpSum(use_vars[j] for j in bpp.BINS), "min_bins"

    for j in bpp.BINS:
      for k in bpp.LIMITS:
        prob += lpSum(bpp.volume[i, k] * assign_vars[i, j] for i in bpp.ITEMS) \
                + waste_vars[j, k] == bpp.capacity[k] * use_vars[j]

    for i in bpp.ITEMS:
        prob += lpSum(assign_vars[i, j] for j in bpp.BINS) == 1

    for i in bpp.ITEMS:
        for j in bpp.BINS:
            prob += assign_vars[i, j] <= use_vars[j]

    for n in range(0, len(bpp.BINS) - 1):
        prob += use_vars[bpp.BINS[n]] >= use_vars[bpp.BINS[n + 1]]

    # Attach the problem data and variable dictionaries to the DipProblem 
    prob.bpp         = bpp
    prob.assign_vars = assign_vars
    prob.use_vars    = use_vars
    prob.waste_vars  = waste_vars

    return prob

def solve(prob):
 
    dippyOpts = {
#                 'doPriceCut' : '1',
                 'CutCGL': '1',
#                'SolveMasterAsIp': '0'
#                'generateInitVars': '1',
#                 'LogDebugLevel': 5,
#                'LogDumpModel': 5,
                 }

    status, message, primals, duals = dippy.Solve(prob, dippyOpts)
  
    if status == LpStatusOptimal:
        return dict((var, var.value()) for var in prob.variables())
    else:
        return None
         
