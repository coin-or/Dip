"""
A formulation in original variables of a wedding seating problem

Authors: Stuart Mitchell 2010
"""
from __future__ import print_function
from builtins import range
import pulp
import argparse

try:
    from src.dippy import DipProblem, DipSolStatOptimal, Solve
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import DipProblem, DipSolStatOptimal, Solve
    from coinor.dippy.examples.gen_func import *

debug_print = False

max_tables = 5
max_table_size = 4
guests = 'A B C D E F G I J K L M N O P Q R'.split()

def happiness(guest_a, guest_b):
    """
    Return the happiness (0 is the best) of allocating two
    guests together in the same table
    """
    return abs(ord(guest_a) - ord(guest_b))

parser = argparse.ArgumentParser(description='Solve a wedding planner problem.')
addDippyArgs(parser)
args = parser.parse_args()

#create the set of possible tables
tables = list(range(max_tables))

possible_seatings = [(g, t) for g in guests
                            for t in tables]

#create a binary variable to model if a guest sits at a particular table
x = pulp.LpVariable.dicts('possible_seatings', possible_seatings, 
                            lowBound = 0,
                            upBound = 1,
                            cat = pulp.LpInteger)

seating_model = DipProblem("Wedding Seating Model (DIP)", pulp.LpMinimize,
                           display_mode = 'off', display_interval = 10000)

#specify the maximum number of guests per table
for table in tables:
    seating_model.relaxation[table] += sum([x[(guest, table)]
                                            for guest in guests]) <= \
                                            max_table_size, \
                                            "Maximum_table_size_%s"%table

#A guest must seated at one and only one table
for guest in guests:
    seating_model += (sum([x[(guest, table)] for table in tables]) == 1, 
                      "Must_seat_%s"%guest)

#create a set of variables to model the objective function
possible_pairs = [(a, b) for a in guests for b in guests if ord(a) < ord(b)]
happy = pulp.LpVariable.dicts('table_happiness', tables,
                              lowBound = 0,
                              upBound = None,
                              cat = pulp.LpContinuous)

seating_model += sum([happy[table] for table in tables])

#create constraints for each possible pair
for table in tables:
    for (a, b) in possible_pairs:
        seating_model.relaxation[table] += \
            happy[table] >= (happiness(a, b) * (x[(a, table)] + 
                                                x[(b, table)] - 1))

def relaxed_solver(prob, table, redCosts, target):
    """
    Generate columns (tables) with negative reduced costs
    """
    dvs = []
    neg_guests = [g for g in guests
                       if redCosts[x[(g,table)]] < 0.0]
    neg_guests.sort()
    # find all possible tables between two end points
    for pos1, pos2 in [(i, j) for i in range(len(neg_guests))
                            for j in range(len(neg_guests))
                                if j > i]:
        # find the suitable guests that can be included in between the end 
        # points
        candidate_guests = [(redCosts[x[(g,table)]], g)
                                    for g in neg_guests[pos1+1:pos2]]
        candidate_guests.sort()
        # pick the best guests (ie those with the negative reduced costs)
        possible_table_inner = [g 
                            for _, g in candidate_guests[:max_table_size-2]]
        #This is the best table between the end points
        possible_table = [neg_guests[pos1]] + possible_table_inner +\
                            [neg_guests[pos2]]
        # calculate the sum of the reduced costs for each of the guests
        neg_cost = sum(redCosts[x[(g, table)]] for g in possible_table)
        table_happiness = happiness(possible_table[0], possible_table[-1])
        rc = neg_cost + table_happiness * redCosts[happy[table]]
        var_values = [(x[(g, table)], 1) 
                      for g in possible_table]
        var_values.append((happy[table], table_happiness))
        dvs.append(dict(var_values))
        if debug_print:
            print('Table: ', table, 'Happiness: ', table_happiness, 'RC: ', rc)
    return DipSolStatOptimal, dvs

#seating_model.relaxed_solver = relaxed_solver

#seating_model.writeLP('wedding_main.lp')
#for table in tables:
#    seating_model.writeRelaxed(table, 'wedding_relax%s.lp' % table);

dippyOpts = addDippyOpts(args)

Solve(seating_model, dippyOpts)

if seating_model.display_mode != 'off':
    numNodes = len(seating_model.Tree.get_node_list())
    if seating_model.Tree.attr['display'] == 'svg':
        seating_model.Tree.write_as_svg(filename = "facility_node%d" % (numNodes + 1), 
                               prevfile = "facility_node%d" % numNodes)
    seating_model.Tree.display()

for table in tables:
    print(table, end=' ')
    for guest in guests:
        if x[(guest,table)].value() >= 0.99:
            print(guest, end=' ')
    print(happy[table].value())

