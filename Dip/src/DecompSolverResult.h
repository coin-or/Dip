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
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//


//===========================================================================//
#ifndef DecompSolverResult_h_
#define DecompSolverResult_h_

//===========================================================================//
/**
 * \class DecompSolverResult
 * \brief Storage of solver result.
 *
 */
//===========================================================================//

//===========================================================================//
#include "Decomp.h"

//===========================================================================//
class DecompSolverResult {

   //----------------------------------------------------------------------//
   /**
    * @name Data.
    * @{
    */
   //----------------------------------------------------------------------//
public:
   int      m_solStatus;
   int      m_solStatus2;
   double   m_objLB;
   double   m_objUB;
   bool     m_isOptimal;
   bool     m_isCutoff;
   int      m_nSolutions;
   double * m_solution;
   /**
    * @}
    */

public:
   /**
    * Default constructors.
    */   
   DecompSolverResult(const int nOrigCols):
      m_solStatus (-1),
      m_solStatus2(-1),
      m_objLB     (-DecompInf),
      m_objUB     ( DecompInf),
      m_isOptimal (false),
      m_isCutoff  (false),
      m_nSolutions(0),
      m_solution  (0)
   {
      m_solution = new double[nOrigCols];
   }
   
   /**
    * Destructor.
    */
   ~DecompSolverResult(){
      UTIL_DELARR(m_solution);
   }
   /**
    * @}
    */
};

#endif
