from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from past.utils import old_div
from .coke_func import *
import importlib as ilib

args = parseArgs()

m = ilib.import_module(args.module)
    
mine_trans = read_table(m.mine_trans_data, int)

cust_trans = read_table(m.cust_trans_data, int, transpose=True)

transport_costs = dict(mine_trans)
transport_costs.update(cust_trans)

cp = CokeProb(supply = m.mine_supply, demand = m.customer_demand,
              LOCATIONS = m.LOCATIONS, build_costs = m.build_costs,
              conversion_factor = m.convert,
              transport_costs = transport_costs)

prob = formulate(cp)

# Set a zero tolerance (Mike Saunders' "magic number")  
prob.tol = pow(pow(2, -24), old_div(2.0, 3.0))

xopt = solve(prob, args)

for l in cp.LOCATIONS:
    for s in cp.SIZES:
        if xopt[prob.buildVars[(l,s)]] > 0:
            print("Build %s %s (%s)" % \
                  (l, s, xopt[prob.buildVars[(l,s)]]))
print()

print_var_table(cp.MINES, cp.LOCATIONS, prob.flowVars)
print()
print_var_table(cp.LOCATIONS, cp.CUSTOMERS, prob.flowVars)

if prob.display_mode != 'off':
    numNodes = len(prob.Tree.get_node_list())
    if ((prob.Tree.attr['display'] == 'pygame') or 
        (prob.Tree.attr['display'] == 'xdot')):
        prob.Tree.display()
