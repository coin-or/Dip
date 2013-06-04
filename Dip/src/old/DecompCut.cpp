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


#include "DecompCut.h"

// --------------------------------------------------------------------- //
bool DecompCut::calcViolation(const CoinPackedVector* row,
                              const double*              x)
{
   //always calculated wrt to original row!
   const double activity = row->dotProduct(x);
   //printf("\nact: %g, m_lb: %g, m_ub: %g",
   //       activity, m_lb, m_ub);
   double violation      = std::max<double>(m_lb - activity, activity - m_ub);
   violation = std::max<double>(0.0, violation);
   //printf("\nviolation = %g", violation);
   setViolation(violation); //should it set it here?
   return violation > 0.0000001;//param?
}

// --------------------------------------------------------------------- //
void DecompCut::print(ostream* os) const
{
   (*os) << "\nCUT"
         << " vio: "    << m_violation
         << " eff: "    << m_effCnt
         << " lb:  "    << m_lb
         << " ub:  "    << m_ub
         << "\n";
   //UtilPrintPackedVector(m_s, os);
   //we don't know anything like in var, we know about s...
   //think....
}

