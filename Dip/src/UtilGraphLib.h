//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

/*
  Author  : Matthew Galati
  Date    : 02/20/03
  Purpose : A class for storing data instances from TSPLIB and VRPLIB.

  02/20/03: Initial version for VRPLIB
  TODO    : TSPLIB
  : apply to KCCP code
*/


#ifndef UTIL_GRAPHLIB_INCLUDED
#define UTIL_GRAPHLIB_INCLUDED

#include "UtilMacros.h"

#include <string>

// ----------------------------------------------------------------------- //
class UtilGraphLib {
private:
   UtilGraphLib(const UtilGraphLib&);
   UtilGraphLib& operator=(const UtilGraphLib&);

public:
   UtilGraphLib() :
      name      (""),
      n_vertices(0),
      n_edges   (0),
      capacity  (0),
      edge_wt   (0),
      vertex_wt (0),
      posx      (0),
      posy      (0),
      coordx    (0),
      coordy    (0),
      coordz    (0) {
   };

   ~UtilGraphLib() {
      UTIL_DELARR(edge_wt);
      UTIL_DELARR(vertex_wt);
      UTIL_DELARR(posx);
      UTIL_DELARR(posy);
      UTIL_DELARR(coordx);
      UTIL_DELARR(coordy);
      UTIL_DELARR(coordz);
   };

public:
   //TODO: make these private and provide access functions?
   //TSPLIB/VRPLIB
   std::string name;
   int      n_vertices;
   int      n_edges;
   int      capacity;
   int*     edge_wt;
   int*     vertex_wt;
   int*     posx;
   int*     posy;
   double* coordx;
   double* coordy;
   double* coordz;

public:
   void read_data(const char* datafile);
   int  compute_icost(const int wtype,
                      const int va,
                      const int vb);
};

#endif
