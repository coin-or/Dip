from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from .bin_pack_func import *
import importlib as ilib

args = parseArgs()

if args.module != None:
    m = ilib.import_module(args.module)
    bpp = BinPackProb(ITEMS = m.ITEMS, volume = m.volume, capacity = m.capacity)
else:
    bpp = BinPackProb(ITEMS  = range(args.numItems),
                      volume = dict(zip(range(args.numItems),[int(i) for i in args.volumes])),
                      capacity = args.capacity)

print(bpp.ITEMS)
print(bpp.volume)
print(bpp.capacity)
  
prob = formulate(bpp, args)

xopt = solve(prob, args)
  
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


