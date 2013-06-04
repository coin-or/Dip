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
// Copyright (C) 2002-2011, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef DECOMP_WAITING_COL_INCLUDE
#define DECOMP_WAITING_COL_INCLUDE

#include "Decomp.h"
#include "DecompVar.h"
#include "UtilMacros.h"
#include "DecompRay.h"



// ---------------------------------------------------------------------- //
class DecompWaitingCol {

 private:
  DecompVar        * m_var;  //s        the variable
  CoinPackedVector * m_col_var;  //(A'' s)  the column
  DecompRay        * m_ray;  // r        the ray
  CoinPackedVector * m_col_ray;  //(A''r) the column


 public:
  inline DecompVar        * getVarPtr() const { return m_var; }
  inline DecompRay        * getRayPtr() const { return m_ray; }
  inline CoinPackedVector * getColVarPtr() const { return m_col_var; }
  inline CoinPackedVector * getColRayPtr() const { return m_col_ray; }
  
  inline const double getVarReducedCost() const  { return m_var->getReducedCost();}
  inline const double getRayReducedCost() const  { return m_ray->getReducedCost();}

  inline const double getVarLowerBound() const   { return m_var->getLowerBound(); }
  inline const double getVarUpperBound() const   { return m_var->getUpperBound(); }

  inline const double getRayLowerBound() const   { return m_ray->getLowerBound(); }
  inline const double getRayUpperBound() const   { return m_ray->getUpperBound(); }

  inline const double getVarOrigCost() const     { return m_var->getOriginalCost(); }
  inline const double getRayOrigCost() const     { return m_ray->getOriginalCost(); }
  
  inline void   deleteColVar() { UTIL_DELPTR(m_col_var); }
  inline void   deleteColRay() { UTIL_DELPTR(m_col_ray); }
  inline void   deleteVar() { UTIL_DELPTR(m_var); }
  inline void   deleteRay() { UTIL_DELPTR(m_ray); }
  inline void   clearVar()  { m_var = 0;          }
  inline void   clearRay()  { m_ray = 0;          }

  inline void   setColVar(CoinPackedVector * col) { m_col_var = col; }
  inline void   setColRay(CoinPackedVector * col) { m_col_ray = col; }

  bool setVarReducedCost(const double      * u, 
			 const DecompStatus    stat);

  bool setRayReducedCost(const double      * u, 
			 const DecompStatus    stat);
  
 public:
  DecompWaitingCol(const DecompWaitingCol & rhs){
    m_var = rhs.m_var;
    m_col_var = rhs.m_col_var;
    m_ray = rhs.m_ray; 
    m_col_ray = rhs.m_col_ray; 
  }
  DecompWaitingCol(DecompVar * var, CoinPackedVector * col) :
    m_var(var),
    m_col_var(col) {}

  DecompWaitingCol(DecompRay * ray, CoinPackedVector * col) :
    m_ray(ray),
    m_col_ray(col) {}
  
  ~DecompWaitingCol() {}
};

#endif
