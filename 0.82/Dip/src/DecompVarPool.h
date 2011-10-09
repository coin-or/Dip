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


#ifndef DECOMP_VAR_POOL_INCLUDE
#define DECOMP_VAR_POOL_INCLUDE

#include "Decomp.h"
#include "DecompWaitingCol.h"

class DecompConstraintSet;

// --------------------------------------------------------------------- //
class is_less_thanD{//member of class instead??
public:
   bool operator()( const DecompWaitingCol & x, 
                    const DecompWaitingCol & y){
      return x.getVarPtr()->getReducedCost() < y.getVarPtr()->getReducedCost();
   }
};

// --------------------------------------------------------------------- //
class DecompVarPool : public std::vector<DecompWaitingCol> {
private:
   DecompVarPool(const DecompVarPool &);
   DecompVarPool & operator=(const DecompVarPool &);

private:
   static const char * classTag;
   bool m_colsAreValid;

public:
   const inline bool colsAreValid() const {return m_colsAreValid;}
   inline void setColsAreValid(bool colsAreValid){
      m_colsAreValid = colsAreValid;
   }

   void print(std::ostream * os = &std::cout) const; //THINK: virtual??
   void reExpand(const DecompConstraintSet & modelCore,
                 const double                tolZero);
   bool isDuplicate(const DecompWaitingCol & wcol);
   bool isDuplicate(const DecompVarList    & vars,
                    const DecompWaitingCol & wcol);
   bool isParallel(const DecompVarList    & vars,
                   const DecompWaitingCol & wcol,
                   const double             maxCosine);
   bool setReducedCosts(const double            * u,
                        const DecompStatus          stat,
                        DecompVarPool::iterator   first,
                        DecompVarPool::iterator   last);
   bool setReducedCosts(const double            * u, 
                        const DecompStatus          stat){
      return setReducedCosts(u, stat, begin(), end());
   }
  
public:
   DecompVarPool() : 
      m_colsAreValid(true) {}

   ~DecompVarPool() {

      //---
      //--- delete any memory that is left in the waiting cols
      //---
      std::vector<DecompWaitingCol>::iterator vi;
      for(vi = begin(); vi != end(); vi++){
         (*vi).deleteVar();
         (*vi).deleteCol();
      }
   }

};

#endif
