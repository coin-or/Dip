#!/usr/bin/env python


__title__   = 'Main module of B&P-cut for the CVPMP and CVPCP'
__version__ = '1.0 Nov 2013'
__author__  = 'Dago Quevedo'
__email__   = 'dago@yalma.fime.uanl.mx'


import  sys
import  draw
import  Global

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
        
import  CVPMP
from    readInstance import *

def main():
    #read instance
    #path=sys.argv[1]
    id,n,p,d,s,w,cxy,V,mytype,id = read('Instances/pmedcap2.dat')
    
    CVPMP.init(n,p,d,s,w,V)
    z, solution = CVPMP.Solver()

    #display solution
    data = Global.path()+'/out.dat'
    f = open(data,'w')
        
    for s in solution:
        i=s[0]
        for j in s[1]:
            if i == j:
                mytype = 1
            else:
                mytype = 0

            f.write('%d\t%d\t%d\t%d\t%d\n'%
                        (mytype, j, i, cxy[j][0], cxy[j][1]))

    f.close();
    draw.draw(data, mytype, id)

main()
