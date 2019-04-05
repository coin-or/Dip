from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from past.utils import old_div
from .coke_func import CokeProb, read_table, formulate, solve, \
                      print_table, print_var_table

if __name__ == '__main__':
    # Python starts here
    convert = 1.3

    mine_supply = {
        "M1":   25.8,
        "M2":  728,
        "M3": 1456,
        "M4":   49,
        "M5":   36.9,
        "M6": 1100,
    }

    LOCATIONS = ["L1", "L2", "L3", "L4", "L5", "L6"]

    build_costs = {
        0:      0,
        75:   4.4,
        150:  7.4,
        225: 10.5,
        300: 13.5,
        375: 16.5,
        450: 19.6,
    }

    customer_demand = {
        "C1":  83,
        "C2":   5.5,
        "C3":   6.975,
        "C4":   5.5,
        "C5": 720.75,
        "C6":   5.5,
    }

    mine_trans_data = """
            L1      L2      L3      L4      L5      L6  
    M1      231737  46813   79337   195845  103445  45186
    M2      179622  267996  117602  200298  128184  49046
    M3      45170   93159   156241  218655  103802  119616
    M4      149925  254305  76423   123534  151784  104081
    M5      152301  205126  24321   66187   195559  88979
    M6      223934  132391  51004   122329  222927  54357
    """

    cust_trans_data = """
            L1      L2      L3      L4      L5      L6 
    C1      6736    42658   70414   45170   184679  111569
    C2      217266  227190  249640  203029  153531  117487
    C3      35936   28768   126316  2498    130317  74034
    C4      73446   52077   108368  75011   49827   62850
    C5      174664  177461  151589  153300  59916   135162
    C6      186302  189099  147026  164938  149836  286307
    """

    mine_trans = read_table(mine_trans_data, int)

    cust_trans = read_table(cust_trans_data, int,
                            transpose=True)

    transport_costs = dict(mine_trans)
    transport_costs.update(cust_trans)
    
    cp = CokeProb(supply = mine_supply, demand = customer_demand,
                  LOCATIONS = LOCATIONS, build_costs = build_costs,
                  conversion_factor = convert,
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

    print_var_table(cp.MINES, LOCATIONS, prob.flowVars)
    print()
    print_var_table(cp.LOCATIONS, cp.CUSTOMERS, prob.flowVars)

    if prob.display_mode != 'off':
        numNodes = len(prob.Tree.get_node_list())
        if ((prob.Tree.attr['display'] == 'pygame') or 
            (prob.Tree.attr['display'] == 'xdot')):
            prob.Tree.display()
