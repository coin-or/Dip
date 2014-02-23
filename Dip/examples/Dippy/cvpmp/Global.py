#!/usr/bin/env python


__title__     = 'General tools'
__version__   = '1.0 Nov 2013'
__author__    = 'Dago Quevedo'
__email__     = 'dago@yalma.fime.uanl.mx'


import math
import os
import sys
import time
#import resource


def euclidean(x1, x2, y1, y2):
    return int(math.sqrt(math.pow(x1 - x2, 2)+math.pow(y1 - y2, 2)))

def indexOf(L, value):
    try:
        i = L.index(value)
    except ValueError:
        i = -1
    
    return i

def path():
    return os.path.abspath(os.path.split(sys.argv[0])[0])