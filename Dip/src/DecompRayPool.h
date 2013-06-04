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


#ifndef DECOMP_RAY_POOL_INCLUDE
#define DECOMP_RAY_POOL_INCLUDE

#include "Decomp.h"
#include "DecompWaitingCol.h"


class DecompConstraintSet;

// --------------------------------------------------------------------- //
/*class is_less_thanD{//member of class instead??
public:
   bool operator()( const DecompWaitingCol & x, 
                    const DecompWaitingCol & y){
      return x.getRayPtr()->getReducedCost() < y.getRayPtr()->getReducedCost();
   }
};
*/
// --------------------------------------------------------------------- //
class DecompRayPool : public std::vector<DecompWaitingCol> {
private:
   DecompRayPool(const DecompRayPool &);
   DecompRayPool & operator=(const DecompRayPool &);

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
   bool isDuplicate(const DecompRayList    & vars,
                    const DecompWaitingCol & wcol);
   bool isParallel(const DecompRayList    & vars,
                   const DecompWaitingCol & wcol,
                   const double             maxCosine);
   bool setReducedCosts(const double            * u,
                        const DecompStatus          stat,
                        DecompRayPool::iterator   first,
                        DecompRayPool::iterator   last);

   bool setReducedCosts(const double            * u, 
                        const DecompStatus          stat){
      return setReducedCosts(u, stat, begin(), end());
   }
  
public:
   DecompRayPool() : 
      m_colsAreValid(true) {}

   ~DecompRayPool() {

      //---
      //--- delete any memory that is left in the waiting cols
      //---
      std::vector<DecompWaitingCol>::iterator vi;
      for(vi = begin(); vi != end(); vi++){
         (*vi).deleteRay();
         (*vi).deleteColRay();
      }
   }

};

#endif
