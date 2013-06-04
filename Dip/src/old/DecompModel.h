//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//


#ifndef DECOMP_MODEL_INCLUDED
#define DECOMP_MODEL_INCLUDED

/*-----------------------------------------------------------------------*/
class DecompConstraintSet;
//class DecompVarList;
//class DecompCutList;
//class DecompVarPool;
//class DecompCutPool;

#include "DecompTypes.h"
#include "DecompVarPool.h"
#include "DecompCutPool.h"

/*-----------------------------------------------------------------------*/
class DecompModel {

private:
   /**
      Disable copy constructors.
   */
   DecompModel(const DecompModel&);
   DecompModel& operator=(const DecompModel&);

public:
   //TODO - change all data members to have m_
   /**
      Model data objects (must be defined by users).
   */
   double*                   objCoeff;      //original c (x-space)
   //DecompConstraintSet     modelCore;     //[A'', b''] : THINK - naming
   //DecompConstraintSet     modelRelax;    //[A',  b' ]

   /**
      Model data objects will be used during algos.
      THINK: belong here or in algos?
   */
   DecompVarList           vars;          //list of vars added to master
   DecompCutList           cuts;
   DecompVarPool           varpool;
   DecompCutPool           cutpool;

public:
   DecompModel() :
      objCoeff(0),
      //modelCore(),
      //modelRelax(),
      vars(),
      cuts(),
      varpool(),
      cutpool()
   {};
   virtual ~DecompModel() {
      UTIL_DELARR(objCoeff);
   }

};

#endif
