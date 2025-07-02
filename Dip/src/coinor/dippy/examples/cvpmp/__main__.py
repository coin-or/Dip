from __future__ import absolute_import
__title__   = 'Main module of B&P-cut for the CVPMP and CVPCP'
__version__ = '1.0 Nov 2013'
__author__  = 'Dago Quevedo'
__email__   = 'dago@yalma.fime.uanl.mx'

from . import  draw
from . import  Global

from .CVPMP import *
from .readInstance import *

def main():

    args = parseArgs()
    
    id,n,p,d,s,w,cxy,V,mytype,id = read(args.instance)
    
    Init(n,p,d,s,w,V)
    z, solution = Solver(args)

    #display solution
    data = 'out.dat'
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
