from __future__ import absolute_import
from builtins import range
#!/usr/bin/env python


__title__   = 'Read instances for the CVPMP and CVPCP'
__version__ = '1.0 Nov 2013'
__author__  = 'Dago Quevedo'
__email__   = 'dago@yalma.fime.uanl.mx'


import string
from . import Global


def read(path):
    file = open(path,'r')
    line = file.readline().split()

    type = int(line[0])
    cxy  = None
    
    #1 - Beasly
    if type == 1:
        id   = int(line[1])
        n    = int(line[2])
        p    = int(line[3])
        line = file.readline().split()
        q    = int(line[0])
        
        cxy  = {}
        d    = {}
        s    = {}
        w    = {}
        V    = [i for i in range(1, n + 1)]
        
        for i in V:
            line    = file.readline().split()
            cxy[i]  = [int(line[1]), int(line[2])]
            s[i]    = q
            w[i]    = int(line[3])
        
        for i in range(1,n+1):
            for j in range(i,n+1):
                e = Global.euclidean(cxy[i][0], cxy[j][0], cxy[i][1], cxy[j][1])
                d[i,j] = e
                d[j,i] = e
    
    
    #2 - GalvaoReVelle
    if type == 2:
        id = int(line[1])
        n  = int(line[2])
        p  = int(line[3])
        
        d  = {}
        s  = {}
        w  = {}
        V  = [i for i in range(1,n+1)]
        
        i = 1
        for s_ in file.readline().split():
            s[i] = int(float(s_))
            i += 1
        
        i = 1
        for w_ in file.readline().split():
            w[i] = int(float(w_))
            i += 1
        
        for i in V:
            j = 1
            for d_ in file.readline().split():
                d[i,j] = int(float(d_))
                j += 1
    
    
    #3 - Lorena
    if type == 3:
        
        id  = int(line[1])
        n   = int(line[2])
        p   = int(line[3])
        
        cxy = {}
        d   = {}
        s   = {}
        w   = {}
        V   = [i for i in range(1,n+1)]
        
        for i in V:
            line    = file.readline().split()
            cxy[i]  = [int(line[0]), int(line[1])]
            s[i]    = int(line[2])
            w[i]    = int(line[3])
        
        for i in range(1,n+1):
            for j in range(i,n+1):
                e = Global.euclidean(cxy[i][0], cxy[j][0], cxy[i][1], cxy[j][1]) 
                d[i,j] = e
                d[j,i] = e

    

    #4 - OR-Library
    if type >= 4 and type <= 7:
        id = int(line[1])
        n  = int(line[2])
        p  = int(line[3])

        
        d  = {}
        s  = {}
        w  = {}
        V  = [i for i in range(1,n+1)]

        for i in V:
            j = 1
            for d_ in file.readline().split():
                d[i,j] = int(float(d_))
                j += 1

        i = 1
        for w_ in file.readline().split():
            w[i] = int(float(w_))
            i += 1
        
        i = 1
        for s_ in file.readline().split():
            s[i] = int(float(s_))
            i += 1
    
    file.close()

    return  id,n,p,d,s,w,cxy,V,id,type
