#!/usr/bin/env python

from mdbin_pack_func import MDBinPackProb, formulate, solve

if __name__ == '__main__':
    # Python starts here
    bpp = MDBinPackProb(ITEMS  = [11, 12, 13, 14, 21, 22],
                        LIMITS = ['CPU', 'RAM'],
                        volume = {(11, 'CPU'): 1,
                                  (12, 'CPU'): 1,
                                  (13, 'CPU'): 1,
                                  (14, 'CPU'): 1,
                                  (21, 'CPU'): 1,
                                  (22, 'CPU'): 1,
                                  (11, 'RAM'): 512,
                                  (12, 'RAM'): 512,
                                  (13, 'RAM'): 512,
                                  (14, 'RAM'): 512,
                                  (21, 'RAM'): 3072,
                                  (22, 'RAM'): 3072,
                                  },
                        capacity = {'CPU': 4, 'RAM': 4096})
  
    prob = formulate(bpp)
  
    prob.tol = pow(pow(2, -24), 2.0 / 3.0)
    xopt = solve(prob)
  
    if xopt is not None:
        for var in prob.variables():
            print var.name, "=", xopt[var]
    else:
        print "Dippy could not find and optimal solution"    
  
    if prob.display_mode != 'off':
        numNodes = len(prob.Tree.get_node_list())
        if ((prob.Tree.attr['display'] == 'pygame') or 
            (prob.Tree.attr['display'] == 'xdot')):
            prob.Tree.display()

  
