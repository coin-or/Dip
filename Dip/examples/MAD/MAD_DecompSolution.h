//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef MAD_DECOMP_SOLUTION_INCLUDED
#define MAD_DECOMP_SOLUTION_INCLUDED

// --------------------------------------------------------------------- //
#include "MAD_DecompApp.h"
#include "DecompSolution.h"

// --------------------------------------------------------------------- //
class MAD_DecompSolution : public DecompSolution {
private:
   const MAD_DecompApp * m_app;

public:
   void print(ostream & os = cout) const {
      int    i, b, border_size;
      double xj;

      const int nOrigRows = m_app->getNOrigRows();
      const int beta      = m_app->getBeta();

      os << "\nBlock Decomposition:";
      vector<unsigned int> border(nOrigRows, 1);
      for(b = 0; b < beta; b++){
         os << "\nBLOCK " << b << ":\t";
         for(i = 0; i < nOrigRows; i++){        
            xj   = m_values[m_app->xIndex(i,b)];
            CoinAssertDebug(UtilIsIntegral(xj));
            CoinAssertDebug(xj <  (1.0 + DecompEpsilon));
            CoinAssertDebug(xj >  (    - DecompEpsilon));
            if(xj > 0.5){
               os << i << " ";
               border[i] = 0;
            }
         }
      }
      border_size = count(border.begin(), border.end(), 1);

      os << "\nBORDER :\t";
      for(i = 0; i < nOrigRows; i++){
         if(!border[i])
            continue;
         os << i << " ";
      }
      os << "\nBORDER Size =  " << border_size << "\n";      
   }
   
private:
   /** @name Copy Constructors
    *
    * Disable the default copy constructors.
    *
    */
   /*
   MAD_DecompSolution(const MAD_DecompSolution &);
   MAD_DecompSolution & operator=(const MAD_DecompSolution &);
   */
 public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes size of solution. */
   MAD_DecompSolution() :
      DecompSolution()
   {
   }
   
   /** Constructor. */
   MAD_DecompSolution(const MAD_DecompApp * app,
                      const int             size,
                      const double        * values,
                      const double          quality) :
      DecompSolution(size, values, quality),
      m_app(app)
   {
   }
   
   virtual ~MAD_DecompSolution() {  
   };
};

#endif
