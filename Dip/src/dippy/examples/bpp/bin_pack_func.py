from builtins import range
from builtins import object
CGL_cuts = False

Bin_antisymmetry = False
Item_antisymmetry = False

Symmetry_branch = False
Most_use_branch = False
Most_assign_branch = False

# Import classes and functions from PuLP
from pulp import LpVariable, lpSum, LpBinary, LpStatusOptimal

# Import any customised paths
try:
    import path
except ImportError:
    pass

# Import dippy (local copy first,
# then a development copy - if python setup.py develop used,
# then the coinor.dippy package
try:
    import src.dippy as dippy
except ImportError:
    import coinor.dippy as dippy

from math import floor, ceil

class BinPackProb(object):
    def __init__(self, ITEMS, volume, capacity):
        self.ITEMS = ITEMS
        self.volume = volume
        self.BINS = list(range(len(ITEMS))) # Create 1 bin for each
                                      # item, indices start at 0
        self.capacity = capacity
    
def formulate(bpp):
    prob = dippy.DipProblem("Bin Packing",
                            display_mode = 'off',
#                           layout = 'bak',
                            display_interval = None,
                            )

    assign_vars = LpVariable.dicts("x",
                                   [(i, j) for i in bpp.BINS
                                    for j in bpp.ITEMS],
                                   cat=LpBinary)
    use_vars    = LpVariable.dicts("y", bpp.BINS, cat=LpBinary)
    waste_vars  = LpVariable.dicts("w", bpp.BINS, 0, None)

    prob += lpSum(waste_vars[i] for i in bpp.BINS), "min_waste"

    for j in bpp.ITEMS:
        prob += lpSum(assign_vars[i, j] for i in bpp.BINS) == 1

    for i in bpp.BINS:
        prob.relaxation[i] += (lpSum(bpp.volume[j] * assign_vars[i, j]
                                for j in bpp.ITEMS) + waste_vars[i] 
                                == bpp.capacity * use_vars[i])

    for i in bpp.BINS:
        for j in bpp.ITEMS:
            prob.relaxation[i] += assign_vars[i, j] <= use_vars[i]

    if Bin_antisymmetry:
        for m in range(0, len(bpp.BINS) - 1):
            prob += use_vars[bpp.BINS[m]] >= use_vars[bpp.BINS[m + 1]]

    if Item_antisymmetry:
        for m in range(0, len(bpp.BINS)):
            for n in range(0, len(bpp.ITEMS)):
                if m > n:
                    i = bpp.BINS[m]
                    j = bpp.ITEMS[n]
                    prob += assign_vars[i, j] == 0

    # Attach the problem data and variable dictionaries
    # to the DipProblem 
    prob.bpp         = bpp
    prob.assign_vars = assign_vars
    prob.use_vars    = use_vars
    prob.waste_vars  = waste_vars

    return prob

def my_branch(prob, sol):

    bounds = None
    
    if Symmetry_branch:
        bounds = symmetry(prob, sol)
  
    if Most_use_branch:
        if bounds is None:
            bounds = most_frac_use(prob, sol)
    
    if Most_assign_branch:
        if bounds is None:
            bounds = most_frac_assign(prob, sol)
    
    return bounds

def my_heuristics(prob, xhat, cost):
#  print "Heuristics..."
    sol = None
  
    if prob.is_root_node:
        prob.is_root_node = False
        if prob.root_heuristic:
            sol = first_fit(prob)
    else:
        if prob.node_heuristic:
            sol = frac_fit(prob, xhat)
  
    if sol is not None:
        return [sol]

def solve(prob, algo = 'PriceCut'):

    if Symmetry_branch or Most_use_branch or Most_assign_branch:
        prob.branch_method = my_branch
#    prob.heuristics = my_heuristics
#    prob.is_root_node = True
#    prob.root_heuristic = True
#    prob.node_heuristic = True
  
    dippyOpts = {}

    if CGL_cuts:
        dippyOpts['CutCGL'] = '1'
    else:
        dippyOpts['CutCGL'] = '0'

    if algo == 'PriceCut':
        dippyOpts['doPriceCut'] = '1'
        dippyOpts['CutCGL'] = '1'
    elif algo == 'Price':
        dippyOpts['doPriceCut'] = '1'
        dippyOpts['CutCGL'] = '0'
    else:
        dippyOpts['doCut'] = '1'

#                'SolveMasterAsIp': '0'
#                'generateInitVars': '1',
#                 'LogDebugLevel': 5,
#                'LogDumpModel': 5,
    dippyOpts['Gurobi'] = {'MipGap':'.05'}

    status, message, primals, duals = dippy.Solve(prob, dippyOpts)
  
    if status == LpStatusOptimal:
        return dict((var, var.value()) for var in prob.variables())
    else:
        return None

def most_frac_use(prob, sol):
    # Get the attached data and variable dicts
    bpp         = prob.bpp
    use_vars    = prob.use_vars
    tol         = prob.tol
  
    most   = float('-inf')
    bin = None
    for j in bpp.BINS:
        alpha = sol[use_vars[j]]
        up   = ceil(alpha)  # Round up to next nearest integer 
        down = floor(alpha) # Round down
        frac = min(up - alpha, alpha - down)
        if frac > tol: # Is fractional?
            if frac > most:
                most = frac
                bin = j
    
    down_lbs = {}
    down_ubs = {}
    up_lbs = {}
    up_ubs = {}
    if bin is not None:
#        print bin, sol[use_vars[bin]]
        down_ubs[use_vars[bin]] = 0.0
        up_lbs[use_vars[bin]] = 1.0
    
        return down_lbs, down_ubs, up_lbs, up_ubs

def most_frac_assign(prob, sol):
    # Get the attached data and variable dicts
    bpp         = prob.bpp
    assign_vars = prob.assign_vars
    tol         = prob.tol
  
    most   = float('-inf')
    assign = None
    for i in bpp.ITEMS:
        for j in bpp.BINS:
            up   = ceil(sol[assign_vars[i, j]])  # Round up to next nearest integer 
            down = floor(sol[assign_vars[i, j]]) # Round down
            frac = min(up - sol[assign_vars[i, j]], sol[assign_vars[i, j]] - down)
            if frac > tol: # Is fractional?
                if frac > most:
                    most = frac
                    assign = (i, j)
    
    down_lbs = {}
    down_ubs = {}
    up_lbs = {}
    up_ubs = {}
    if assign is not None:
#    print assign, sol[assign_vars[assign]]
        down_ubs[assign_vars[assign]] = 0.0
        up_lbs[assign_vars[assign]] = 1.0
    
        return down_lbs, down_ubs, up_lbs, up_ubs

def symmetry(prob, sol):
    # Get the attached data and variable dicts
    bpp      = prob.bpp
    use_vars = prob.use_vars
    tol      = prob.tol
  
    alpha = sum(sol[use_vars[j]] for j in bpp.BINS)
#      print "# bins =", alpha
    up    = int(ceil(alpha))  # Round up to next nearest integer 
    down  = int(floor(alpha)) # Round down
    frac  = min(up - alpha, alpha - down)
    if frac > tol: # Is fractional?
#    print "Symmetry branch"
    
        down_lbs = {}
        down_ubs = {}
        up_lbs = {}
        up_ubs = {}
        for n in range(up - 1, len(bpp.BINS)):
            down_ubs[use_vars[bpp.BINS[n]]] = 0.0
#           print down_ubs
        for n in range(up): # Same as range(0, up)
            up_lbs[use_vars[bpp.BINS[n]]] = 1.0
#           print up_lbs

        return down_lbs, down_ubs, up_lbs, up_ubs
  
def fit(prob, order):
    bpp         = prob.bpp
    use_vars    = prob.use_vars
    assign_vars = prob.assign_vars
    waste_vars  = prob.waste_vars
    tol         = prob.tol

    sol = {}

    for j in bpp.BINS:
        sol[use_vars[j]] = 1.0
        for i in bpp.ITEMS:
            sol[assign_vars[i, j]] = 0.0
        sol[waste_vars[j]] = bpp.capacity

    assigned = {}
    for i in bpp.ITEMS:
        assigned[i] = False

    for (i, j) in order:
        if (not assigned[i]) and (bpp.volume[i] <= sol[waste_vars[j]]):
            assigned[i] = True
            sol[assign_vars[i, j]] = 1.0
            sol[waste_vars[j]] -= bpp.volume[i]
  
    for j in bpp.BINS:
        if sol[waste_vars[j]] > bpp.capacity - tol:
            sol[use_vars[j]] = 0.0
            sol[waste_vars[j]] = 0.0
      
#   print sol
  
    return sol

import operator

def first_fit(prob):
#  print "first fit..."
  
    bpp = prob.bpp
    sorted_volume = sorted(iter(bpp.volume.items()), key=operator.itemgetter(1), reverse=True)
    sorted_ITEMS = [i for (i, v) in sorted_volume]
  
#   print sorted_ITEMS

    order = [(i, j) for i in sorted_ITEMS for j in bpp.BINS]
  
#   print order
  
    sol = fit(prob, order)
#   print sol
    return sol

def frac_fit(prob, xhat):
#    print "frac fit..."
  
    bpp = prob.bpp
    assign_vars = prob.assign_vars
  
    assign = dict(((i, j), xhat[assign_vars[i, j]]) for i in bpp.ITEMS for j in bpp.BINS)
#   print assign
  
    sorted_assign = sorted(iter(assign.items()), key=operator.itemgetter(1), reverse=True)
    order = [(i, j) for ((i, j), x) in sorted_assign]
   
#   print order
  
    sol = fit(prob, order)
#   print sol
    return sol
         
