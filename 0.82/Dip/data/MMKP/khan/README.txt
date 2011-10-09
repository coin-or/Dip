This is the format for benchmark data file for the MMKP. The comments are inside the parantheses. The index of an item in a group are numbered from 0 to (l-1). Please see the original paper for the symbols used in the format.

Solving the Knapdack Problem for Adaptive Multimedia System, Studia Informatica Universalis, An International Journal on Informatics, Special Issue on Cutting, Packing and Knapsacking Problems, vol. 2, No 1, pp. 161-182, 2002.




Beginning of the files (denoted MMKP1, ..., MMKP13) from the next line:
n  l m 
 R_1  R_2  ......... R_m (Available resources: from resource 1 to resource m)
 1	(Data of group 1)
 v_11 (value of item 1)  r_111 r_112  .....  r_11m (Resource consumption by item 1)
 v_12 (value of item 2)  r_211 r_212  .....  r_21m (Resource consumption by item 2)
 .
 .
 v_1l (value of item l)  r_1l1 r_1l2  .....  r_1lm (Resource consumption by item l)
 2
 .........................
 .........................
 .........................
 n
 .........................
 .........................
 ......................... 


--------------------------------------------------------------------------------------
 Solutions by  Mosers Solution :
========================
 \rho_1  \rho_2  ......... \rho_l (index of groups)  V_MOSER(Solution value by Moser)


--------------------------------------------------------------------------------------
 Solutions by  HEU:
========================
 \rho_1  \rho_2  ......... \rho_l (index of groups)  V_HEU(Solution value by HEU)
 Upper bound


--------------------------------------------------------------------------------------
 V_UPPER(upper bound of the problem):
===============================
 Solutions by  Exact Solution:
===============================
 \rho_1  \rho_2  ......... \rho_l (index of groups)  V_UPPER(Solution value by UPPER)
--------------------------------------------------------------------------------------
