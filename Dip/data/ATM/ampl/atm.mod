set A;
set D;
set AD within {A,D};

param Ca{AD};
param Cb{AD};
param Cc{AD};
param Cd{AD};
param Ce{AD};
param Cw{AD};

param K{A};
param B{D};

var x1{A}  >= 0 <= 1;
var x2{A}  >= 0 <= 1;
var x3{A}  >= 0;
var fp{AD} >= 0;
var fm{AD} >= 0;
var v{AD}  binary;

minimize cost: 
   sum {(a,d) in AD} (fp[a,d] + fm[a,d]);

budget{d in D}:
   sum {(a,d) in AD} (fp[a,d] - fm[a,d]) <= B[d];

fdefine{(a,d) in AD}: 
   fp[a,d] - fm[a,d] =
     Ca[a,d] * x1[a]       + 
     Cb[a,d] * x2[a]       + 
     Cc[a,d] * x1[a]*x2[a] +
     Cd[a,d] * x3[a]       + 
     Ce[a,d];

link{(a,d) in AD}:
  fm[a,d] <= Cw[a,d] * v[a,d];

count{a in A}:
   sum {(a,d) in AD} v[a,d] <= K[a]; 
