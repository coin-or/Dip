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

/*===========================================================================*
 * This file is part of the Abstract Library for Parallel Search (ALPS).     *
 *                                                                           *
 * ALPS is distributed under the Common Public License as part of the        *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors: Yan Xu, SAS Institute Inc.                                       *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2006, Lehigh University, Yan Xu, and Ted Ralphs.       *
 *===========================================================================*/

#include <iomanip>
#include <iostream>
#include <set>


#include "BcpsDecompModel.h"
#include "BcpsDecompSolution.h"

//#############################################################################

void 
BcpsDecompSolution::print(std::ostream& os) const
{
    os <<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);

    os << "-------------------------" <<std::endl;
    for (int i = 0; i < size_; ++i) {
	if (fabs(value_[i]) > 1.0e-7) {
	    os << std::setw(6) << i << " " << value_[i] << std::endl;
	}
    }
    os << "-------------------------" <<std::endl;
    os <<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);

}

//can't this be inherited?
/** The method that encodes the node into a encoded object. */
#if 0
AlpsEncoded*
BcpsDecompSolution::encode() const 
{ 
  //  BcpsDecompEncoded* encoded = new BcpsDecompEncoded(typeid(*this).name());
  AlpsEncoded* encoded = new AlpsEncoded("ALPS_SOLUTION");
  
  encoded->writeRep(size_);     // Base operand of `->' has non-pointer type
  encoded->writeRep(value_, size_);
  encoded->writeRep(objective_);

  return encoded; 
}

//#############################################################################

/** The method that decodes the node from a encoded object. */
// Note: write and read sequence MUST same! 
AlpsKnowledge* 
BcpsDecompSolution::decode(AlpsEncoded& encoded) const 
{ 
  int s;
  double obj;
  double* val = 0;
  encoded.readRep(s);        
  encoded.readRep(val, s);        // s must immediately before sol
  encoded.readRep(obj);
  BcpsDecompSolution* sol = new BcpsDecompSolution(s, val, obj);

  if (val != 0) {
      delete [] val;
      val = 0;
  }
  
  return sol;
}
#endif
//#############################################################################
