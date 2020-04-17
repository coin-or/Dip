from __future__ import absolute_import
__title__   = 'Main module of B&P-cut for the CVPMP and CVPCP'
__version__ = '1.0 Nov 2013'
__author__  = 'Dago Quevedo'
__email__   = 'dago@yalma.fime.uanl.mx'

import os, sys
from os.path import join, dirname
from inspect import getfile
import coinor.dippy.examples.cvpmp

from . import  draw
from . import  Global

from . import  CVPMP
from .readInstance import *

def main():
    if len(sys.argv) > 1:
        if sys.argv[1] == '-h' or sys.argv[1] == '--help' or len(sys.argv) > 2:
            print('Usage: coke <module_name>')
            print('       module_name : Python module containing instance data')
            print('                     For example file, check directory')
            print('                    ', dirname(getfile(coinor.dippy.examples.coke)))
            exit()
        else:
            path = sys.argv[1]
    else:
        path = join(dirname(getfile(coinor.dippy.examples.cvpmp)), 'Instances','pmedcap1.dat')

    id,n,p,d,s,w,cxy,V,mytype,id = read(path)
    
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
