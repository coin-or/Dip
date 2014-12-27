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

#ifndef AlpsDecompSolution_h
#define AlpsDecompSolution_h

//===========================================================================//
#include "AlpsSolution.h"
#include "AlpsDecompModel.h"

//===========================================================================//
class AlpsDecompSolution : public AlpsSolution {
protected:
   /** Length of solution (number of columns). */
   int      m_size;

   /** Solution values. */
   double* m_values;

   /** Quality of solution (bound wrt to objective). */
   double   m_quality;

   /** Pointer to DecompApp for the print function. */
   const DecompApp* m_app;

public:
   /** @name Helper functions (public). */

   /** Get length of solution. */
   inline const int getSize() const {
      return m_size;
   }

   /** Get solution values. */
   inline const double* getValues() const {
      return m_values;
   }

   /** Get quality of solution. */
   inline const double getQuality() const {
      return m_quality;
   }

   const DecompApp* getModel() const {return m_app;} 

public:
   AlpsDecompSolution(const DecompApp* m) :
      AlpsSolution(),
      m_size      (0),
      m_values    (0),
      m_quality   (1e75),
      m_app       (m) {}


   AlpsDecompSolution(const int             size,
                      const double*         values,
                      const double          quality,
                      const DecompApp*      app = NULL,
                      const int             depth = -1,
                      const AlpsNodeIndex_t index = -1) :
      AlpsSolution(index, depth),
      m_size      (size),
      m_values    (0),
      m_quality   (quality),
      m_app       (app) {
      CoinAssert(m_size > 0);
      m_values = new double[m_size];
      CoinAssertHint(m_values, "Error: Out of Memory");
      memcpy(m_values, values, sizeof(double) * m_size);
   }

   virtual ~AlpsDecompSolution() {
      UTIL_DELARR(m_values);
   };

   /** Print out the solution.*/
   virtual void print(std::ostream& os) const {
      if (m_app) {
         DecompAlgo*           decompAlgo = m_app->getDecompAlgo();
         DecompConstraintSet* modelCore
         = decompAlgo->getModelCore().getModel();
         m_app->printOriginalSolution(m_size,
                                      modelCore->getColNames(),
                                      m_values);
      }

      /*int i;
      os << setiosflags(ios::fixed|ios::showpoint)
         << setw(14);

      os << "-------------------------" << endl;
      os << "Quality = " << getQuality() << endl;
      for (i = 0; i < m_size; i++) {
         if (!UtilIsZero(m_values[i])){
       os << setw(6) << i << " " << m_values[i] << endl;
         }
      }
      os << "-------------------------" << endl;
      os << resetiosflags(ios::fixed|ios::showpoint|ios::scientific); */
   }


   /** The method that encodes the node into an encoded object **/
   virtual AlpsEncoded* encode() const {
      AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeSolution);
      encoded->writeRep(m_size);
      encoded->writeRep(m_values, m_size);
      encoded->writeRep(m_quality);
      return encoded;
   }

   /** The method that decodes the node from an encoded object **/
   virtual AlpsKnowledge* decode(AlpsEncoded& encoded) const {
      int s, q;
      double* v = 0;
      encoded.readRep(s);
      encoded.readRep(v, s);
      encoded.readRep(q);
      //TODO: change of design of AlpsDecompSolution
      //       constructor
      return new AlpsDecompSolution(s, v, q, getModel());
   }


};

#endif