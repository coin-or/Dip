#!/usr/bin/env python

from __future__ import division
from builtins import range
from past.utils import old_div
__title__     = 'B&P-cut for the Capacitated Vertex p-Median Problem (CVPMP)'
__version__   = '1.0 Nov 2013'
__author__    = 'Dago Quevedo'
__email__     = 'dago@yalma.fime.uanl.mx'

import  sys
from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize
from    math import *
try:
    import path
except ImportError:
    pass
        
try:
    import src.dippy as dippy
    from src.dippy import DipSolStatOptimal
except ImportError:
    import coinor.dippy as dippy
    from coinor.dippy import DipSolStatOptimal

#Globar vars

n  = None
p  = None
d  = None
s  = None
w  = None
V  = None
x  = None
y  = None

tol          = pow(pow(2, -24), old_div(2.0, 3.0))
display_mode = 'off'

def init(_n,_p,_d,_s,_w,_V):
    
    global n, p, d, s, w, V
    n   = _n
    p   = _p
    d   = _d
    s   = _s
    w   = _w
    V   = _V

def solve_subproblem(prob, i, redCosts, target):
    
    vars   = [x[(i, j)] for j in V]
    obj    = [max(-redCosts[x[(i, j)]], 0) for j in V]
    weights= [w[j] for j in V]

    #Solver a knapsack for i
    z, solution = KP01(obj, weights, s[i])
    rc = redCosts[y[i]] - z
    
    #Cost
    cost           = sum([d[i,V[j]] for j in solution])

    #Cluster of customers for i
    var_val        = dict([(vars[j], 1) for j in solution])
    var_val[y[i]]  = 1

    return DipSolStatOptimal, [var_val]

def KP01(obj, weights, capacity):
    
    assert len(obj) == len(weights)

    n = len(obj)
    if n == 0:
        return 0, []
    
    c     = [[0]*(capacity+1) for i in range(n)]
    added = [[False]*(capacity+1) for i in range(n)]
    
    for i in range(n):
        for j in range(capacity + 1):
            if (weights[i] > j):
                c[i][j] = c[i-1][j]
            else:
                c_add = obj[i] + c[i-1][j-weights[i]]
                if c_add > c[i-1][j]:
                    c[i][j]     = c_add
                    added[i][j] = True
                else:
                    c[i][j] = c[i-1][j]
    
    i = n-1
    j = capacity
    
    solution = []
    while i >= 0 and j >= 0:
        if added[i][j]:
            solution.append(i)
            j -= weights[i]
        i -= 1
    
    return c[n-1][capacity], solution



def Solver():
    
    global  x,y

    prob = dippy.DipProblem("CVPMP", display_mode = display_mode,
                            layout = 'dot', display_interval = 0)
    
    X = [(i, j) for i in V for j in V]
    x = LpVariable.dicts("x", X, 0, 1, LpBinary)
    y = LpVariable.dicts("y", V, 0, 1, LpBinary)

    prob += (lpSum(d[i, j] * x[(i, j)] for i in V for j in V), "min")
    
    #linking constraints
    for j in V:
        prob += lpSum(x[(i,j)] for i in V) == 1
    
    #non-relaxing
    for i in V:
        prob.relaxation[i] += lpSum(w[j] * x[(i, j)]
                                    for j in V) <= s[i]*y[i]

    prob += lpSum(y[i] for i in V) == p

    prob.relaxed_solver = solve_subproblem

    dippy.Solve(prob, {
        'TolZero'           : '%s' % tol,
        'doCut'        : '1',
        'generateInitVars'  : '1',
        'CutCGL'            : '1',
    })


    #Make solution
    solution = []
    for i in V:
        if y[i].varValue:
            cluster = []
            for j in V:
                if x[(i, j)].varValue:
                    cluster.append(j)

            solution.append((i,cluster))


    return round(prob.objective.value()), solution

