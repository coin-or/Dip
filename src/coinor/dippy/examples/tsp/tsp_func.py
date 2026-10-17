from __future__ import print_function
from builtins import range
from pulp import *
import argparse
import importlib as ilib

try:
    from src.dippy import DipProblem, DipSolStatOptimal
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import DipProblem, DipSolStatOptimal
    from coinor.dippy.examples.gen_func import *

def parseArgs():
    
    parser = argparse.ArgumentParser(
        description='Solve a traveling salesman problem.')
    parser.add_argument('--module', '-m', metavar = 'module name', 
                        help='name of the Python module from which to import data',
                        default = 'coinor.dippy.examples.tsp.tsp_data')

    addDippyArgs(parser)

    args = parser.parse_args()

    return(args)
 
def formulate(module_name):

    m = ilib.import_module(module_name)
    
    prob = DipProblem()

    arc_vars = LpVariable.dicts("UseArc", m.ARCS, 0, 1, LpBinary)

    # objective
    prob += lpSum(m.ARC_COSTS[x] * arc_vars[x] for x in m.ARCS)

    # degree constraints
    for city in m.CITIES:
        prob += lpSum(arc_vars[x] for x in m.CITY_ARCS[city]) \
                == 2

    # dictionary for symmetric arcs, can be
    # accessed using (i, j) and (j, i)
    symmetric = {}
    for i in m.CITIES:
        for j in m.CITIES:
            if i < j:
                symmetric[(i, j)] = (i, j)
                symmetric[(j, i)] = (i, j)

    prob.CITIES = m.CITIES
    prob.symmetric = symmetric
    prob.arc_vars = arc_vars

    return prob

# generate subtour elimination constraints
def generate_cuts(prob, sol):
    cons = []
    not_connected = set(prob.CITIES)

    while not_connected:
        start = not_connected.pop()
        nodes, arcs = get_subtour(prob, sol, start)
        if len(nodes) == len(arcs) and \
           len(nodes) < len(prob.CITIES):
            cons.append( lpSum(prob.arc_vars[a] for a in arcs) \
                         <= len(arcs) - 1 )
        nodes.remove(start)
        not_connected -= nodes

    return cons

def is_solution_feasible(prob, sol, tol):
    nodes, arcs = get_subtour(prob, sol, 0)

    return len(nodes) == len(arcs) and \
           len(nodes) == len(prob.CITIES)

def get_subtour(prob, sol, node):
    # returns: list of nodes and arcs
    # in subtour containing node
    nodes = set([node])
    arcs = set([])
    not_processed = set(prob.CITIES)
    to_process = set([node])

    tol = 1e-6
    one = 1 - tol

    while to_process:
        c = to_process.pop()
        not_processed.remove(c)
        new_arcs = [ prob.symmetric[(c, i)] \
                     for i in not_processed \
                     if sol[ \
                         prob.arc_vars[prob.symmetric[(c, i)]]]
                     > one]                     
        new_nodes = [ i for i in not_processed \
                      if prob.symmetric[(i, c)] in new_arcs ]
        arcs |= set(new_arcs)
        nodes |= set(new_nodes)
        to_process |= set(new_nodes)

    return nodes, arcs

