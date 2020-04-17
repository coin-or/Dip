from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from past.utils import old_div
from .coke_func import *
import sys
import importlib as ilib
from os.path import dirname
from inspect import getfile
import coinor.dippy.examples.coke

if len(sys.argv) > 1:
    if sys.argv[1] == '-h' or sys.argv[1] == '--help' or len(sys.argv) > 2:
        print('Usage: coke <module_name>')
        print('       module_name : Python module containing instance data')
        print('                     For example file, check directory')
        print('                    ', dirname(getfile(coinor.dippy.examples.coke)))
        exit()
    else:
        module_name = sys.argv[1]
else:
    module_name = 'coinor.dippy.examples.coke.coke_data'

m = ilib.import_module(module_name)
    
mine_trans = read_table(m.mine_trans_data, int)

cust_trans = read_table(m.cust_trans_data, int,
                        transpose=True)

transport_costs = dict(mine_trans)
transport_costs.update(cust_trans)

cp = CokeProb(supply = m.mine_supply, demand = m.customer_demand,
              LOCATIONS = m.LOCATIONS, build_costs = m.build_costs,
              conversion_factor = m.convert,
              transport_costs = transport_costs)

prob = formulate(cp)

# Set a zero tolerance (Mike Saunders' "magic number")  
prob.tol = pow(pow(2, -24), old_div(2.0, 3.0))

xopt = solve(prob)

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
