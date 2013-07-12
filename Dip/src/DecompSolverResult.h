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
// Copyright (C) 2002-2011, Lehigh University, Matthew Galati, and Ted Ralphs//
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
#include "DecompSolution.h"

//===========================================================================//
class DecompSolverResult {

   //----------------------------------------------------------------------//
   /**
    * @name Data.
    * @{
    */
   //----------------------------------------------------------------------//
public :
   int       m_solStatus;
   int       m_solStatus2;
   double    m_objLB;
   double    m_objUB;
   bool      m_isOptimal;
   bool      m_isCutoff;
   
   int       m_nPoints;
   int       m_nRays; 

   bool      m_isUnbounded; 

   std::vector< std::vector<double> > m_point;

   std::vector< std::vector<double> > m_ray; 
   /**
    * @}
    */

public:
   const double * getSolution(const int solIndex){
      std::vector<double> & solution = m_point[solIndex];
      return &solution[0];
   }

   /**
    * Default constructors.
    */   
   DecompSolverResult():
      m_solStatus (-1),
      m_solStatus2(-1),
      m_objLB     (-DecompInf),
      m_objUB     ( DecompInf),
      m_isOptimal (false),
      m_isCutoff  (false),	
      m_nPoints(0),
      m_nRays(0)
   {
   }

   DecompSolverResult(const DecompSolution * solution):
      m_solStatus (-1),
      m_solStatus2(-1),
      m_objLB     (-DecompInf),
      m_objUB     ( DecompInf),
      m_isOptimal (false),
      m_isCutoff  (false),
      m_nPoints(0),
      m_nRays(0)	
   {
      const double * values = solution->getValues();
      if(!values)
	 return;
      m_nPoints = 1;
      m_nRays = 0; 

      m_objUB      = solution->getQuality();
      std::vector<double> sol(values, values + solution->getSize());
      m_point.push_back(sol);
   }
   
   /**
    * Destructor.
    */
   ~DecompSolverResult(){
   }
   /**
    * @}
    */
};

#endif
