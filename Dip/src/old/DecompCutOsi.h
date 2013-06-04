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


#ifndef DECOMP_CUTOSI_HPP
#define DECOMP_CUTOSI_HPP


#include "UtilHash.h"
#include "DecompCut.h"
#include "OsiRowCut.hpp"

//THINK:
//why?? if we really want to be able to let the user use OsiRowCut's
//then just make this public OsiRowCut or something... then need multiple
//inheritance? how does ABC do?

//really want a generic cut implementation here...
//or just use osi?


class DecompCutOsi : public DecompCut {
private:
   DecompCutOsi(const DecompVar&);
   DecompCutOsi& operator=(const DecompVar&);

private:
   OsiRowCut m_osiCut;
   /* THINK: seems a waste to have to copy construct the OsiRowCut, a pointer
      should be enough - but error on pure virtual */
   /* at least make it a pointer to OsiRowCut? */
public:
   //is this an expensive operation?
   //temp fix
   char sense() const {
      double lb_ = m_osiCut.lb();
      double ub_ = m_osiCut.ub();

      if      ( lb_ == ub_ ) {
         return 'E';
      } else if ( lb_ == -DecompInf && ub_ == DecompInf ) {
         return 'N';
      } else if ( lb_ == -DecompInf ) {
         return 'L';
      } else if ( ub_ == DecompInf ) {
         return 'G';
      } else {
         return 'R';
      }
   }

   double rhs() const {
      double lb_ = m_osiCut.lb();
      double ub_ = m_osiCut.ub();

      if      ( lb_ == ub_ ) {
         return ub_;
      } else if ( lb_ == -DecompInf && ub_ == DecompInf ) {
         return 0.0;
      } else if ( lb_ == -DecompInf ) {
         return ub_;
      } else if ( ub_ == DecompInf ) {
         return lb_;
      } else {
         return ub_;
      }
   }

   void setStringHash() {
      //we cannot trust osi row cuts sense, since cpx and clp have different infinities...
      m_strHash = UtilCreateStringHash(m_osiCut.row().getNumElements(),
                                       m_osiCut.row().getIndices(),
                                       m_osiCut.row().getElements(),
                                       //m_osiCut.sense(),
                                       sense(),
                                       //m_osiCut.rhs()
                                       rhs()
                                      );
      //ranges?
   }
   void setStringHash(CoinPackedVector* row) {
      m_strHash = UtilCreateStringHash(row->getNumElements(),
                                       row->getIndices(),
                                       row->getElements(),
                                       //m_osiCut.sense(),
                                       sense(),
                                       //m_osiCut.rhs()
                                       rhs()
                                      );
      //ranges?
   }

   void setBounds() {
      setLowerBound(m_osiCut.lb());
      setUpperBound(m_osiCut.ub());
   }

   //think about when is this used?
   void expandCutToRow(CoinPackedVector* row) {
      row->setVector(m_osiCut.row().getNumElements(),
                     m_osiCut.row().getIndices(),
                     m_osiCut.row().getElements(),
                     DECOMP_TEST_DUPINDEX);
      /* TODO: tests for dups by default - shut this off for production */
   }

public:
   void print(ostream* os = &cout) const {
      (*os).precision(2);
      (*os) << endl;
      const int*     ind = m_osiCut.row().getIndices();
      const double* els = m_osiCut.row().getElements();

      for (int i = 0; i < m_osiCut.row().getNumElements(); i++) {
         (*os) << " + " << els[i] << " x[" << ind[i] << "]";
      }

      if (getLowerBound() < -1.0e10 / 2) { //INF?
         (*os) << " lb: -INF";
      } else {
         (*os) << " lb: " << getLowerBound();
      }

      if (getUpperBound() > 1.0e10 / 2) { //INF?
         (*os) << " ub: INF";
      } else {
         (*os) << " ub: " << getUpperBound();
      }

      (*os) << " vio: " << getViolation() << "\n";
   }

public:
   DecompCutOsi(OsiRowCut& osiCut)
      : DecompCut(), m_osiCut(osiCut) {
      setBounds();
   };
   virtual ~DecompCutOsi() {}

};

#endif
