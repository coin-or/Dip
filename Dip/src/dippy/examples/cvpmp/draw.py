#!/usr/bin/env python


from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from past.utils import old_div
__title__   = 'Display graphic solution for the CVPMP and CVPCP'
__version__ = '1.0 Nov 2013'
__author__  = 'Dago Quevedo'
__email__   = 'dago@yalma.fime.uanl.mx'


import  string
from    tkinter import  * 

def draw(path, mytype, id):
    f = open(path,"r")
    V = []
    P = []
    
    for line in f:
        Z = line.split()
        V.append([int(Z[0]),int(Z[1]),int(Z[2]),float(Z[3]),float(Z[4])])
        if int(Z[0]) ==  1:
            P.append([int(Z[0]),int(Z[1]),int(Z[2]),float(Z[3]),float(Z[4])])
    
    f.close()

    x_max = max(V,key = lambda x:x[3])[3]
    y_max = max(V,key = lambda x:x[4])[4]
    x_min = min(V,key = lambda x:x[3])[3]
    y_min = min(V,key = lambda x:x[4])[4]
    
    #Normalizacion de valores

    if mytype == 1:
        scale = 600
        delta = 20

    if mytype == 0:
        scale = 650
        delta = 20
    

    for i in range(len(V)):
        V[i][3] = ((old_div((V[i][3] - x_min), (x_max - x_min))) * scale) + delta
        V[i][4] = ((old_div((V[i][4] - y_min), (y_max - y_min))) * scale) + delta
        if i  <  len(P):
            P[i][3] = ((old_div((P[i][3] - x_min), (x_max - x_min))) * scale) + delta
            P[i][4] = ((old_div((P[i][4] - y_min), (y_max - y_min))) * scale) + delta

    _d_ = 0
    
    C = ["#87CEFA","#C0C0C0","#FFA500","#DDA0DD","#9ACD32", "#9ACD32",
         "#E99699","#B8F7B8","#D8D8F9","#FAD6A5","#F4FBA6", "#FF434F",
         "#87CEFA","#C0C0C0","#FFA500","#DDA0DD","#9ACD32", "#9ACD32",
         "#E99699","#B8F7B8","#D8D8F9","#FAD6A5","#F4FBA6", "#FF434F",
         "#87CEFA","#C0C0C0","#FFA500","#DDA0DD","#9ACD32", "#9ACD32",
         "#E99699","#B8F7B8","#D8D8F9","#FAD6A5","#F4FBA6", "#FF434F",
         "#87CEFA","#C0C0C0","#FFA500","#DDA0DD","#9ACD32", "#9ACD32",
         "#E99699","#B8F7B8","#D8D8F9","#FAD6A5","#F4FBA6", "#FF434F",
         "#87CEFA","#C0C0C0","#FFA500","#DDA0DD","#9ACD32", "#9ACD32",
         "#E99699","#B8F7B8","#D8D8F9","#FAD6A5","#F4FBA6", "#FF434F",
         "#87CEFA","#C0C0C0","#FFA500","#DDA0DD","#9ACD32", "#9ACD32",
         "#E99699","#B8F7B8","#D8D8F9","#FAD6A5","#F4FBA6", "#FF434F"]
    
    root = Tk()
    root.wm_attributes("-topmost", 1)

    root.title("Solution " + str(type) + " - " + str(id))
    canvas = Canvas(root,width = scale + (delta * 2), height = scale + (delta * 2), 
                    bg = 'white')  
    canvas.pack(expand = YES, fill = BOTH)

    for _k_ in range(len(P)):
        k = P[_k_]
        for j in V:
            if j[2] == k[1]:
                if j[0] == 0:
                    canvas.create_line(
                        k[3] + 10,k[4] + 10, 
                        j[3] + 7, j[4] + 7, 
                        dash = (3,3))
                else:
                    canvas.create_line(
                        k[3] + 10,k[4] + 10, 
                        j[3] + 7,j[4] + 7,fill = "red",width = 1.5, 
                        dash = (3,3))
                

                if j[1] > 9:_d_ = 2
                else     :_d_ = 0
                canvas.create_oval(
                        j[3],j[4], 
                        j[3] + 15,j[4] + 15, 
                                   width = 1, fill = C[_k_])
                
                canvas.create_text(
                        j[3] + (5 - _d_),
                        j[4] + 7, 
                                   text = str(j[1]),anchor = "w",fill = "black", 
                                   font = ("Arial", 10))
        canvas.create_rectangle(
                k[3],k[4], 
                k[3] + 20,k[4] + 20, 
                width = 2, fill = C[_k_]) 

        if k[1] > 9:_d_ = 3
        else       :_d_ = 0
        canvas.create_text(
                 k[3] + (6 - _d_),k[4] + 10, 
                 text = str(k[1]),anchor = "w",fill = "black", 
                 font = ("Arial", 14))

    root.mainloop()

