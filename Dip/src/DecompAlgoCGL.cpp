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
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "DecompAlgo.h"
#include "DecompCutOsi.h"
//===========================================================================//
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CoinPackedMatrix.hpp"

using namespace std;

//===========================================================================//
int DecompAlgoCGL::initGenerators(const int doClique,
                                  const int doOddHole,
                                  const int doFlowCover,
                                  const int doKnapCover,
                                  const int doMixIntRound,
                                  const int doGomory)
{
   int status = DecompStatOk;

   if (doClique) {
      m_genClique = new CglClique;
      m_genClique->setStarCliqueReport(false);
      m_genClique->setRowCliqueReport (false);

      if (!m_genClique) {
         return DecompStatOutOfMemory;
      }
   }

   if (doOddHole) {
      m_genOddHole = new CglOddHole;

      if (!m_genClique) {
         return DecompStatOutOfMemory;
      }
   }

   if (doFlowCover) {
      m_genFlowCover = new CglFlowCover;

      if (!m_genFlowCover) {
         return DecompStatOutOfMemory;
      }
   }

   if (doKnapCover) {
      m_genKnapCover = new CglKnapsackCover;

      if (!m_genKnapCover) {
         return DecompStatOutOfMemory;
      }
   }

   if (doMixIntRound) {
      m_genMixIntRound = new CglMixedIntegerRounding2;

      if (!m_genMixIntRound) {
         return DecompStatOutOfMemory;
      }
   }

   if (doGomory) {
      m_genGomory = new CglGomory;

      if (!m_genGomory) {
         return DecompStatOutOfMemory;
      }
   }

   return status;
}

//===========================================================================//
int DecompAlgoCGL::generateCuts(OsiSolverInterface* cutGenSI,
                                OsiSolverInterface* masterSI,
                                double*              xhat,
                                vector<int>&         integerVars,
                                DecompCutList&       newCuts)
{
   OsiCuts      osiCuts;
   int          status           = DecompStatOk;
   int          nCliqueCuts      = 0;
   int          nOddHoleCuts     = 0;
   int          nFlowCoverCuts   = 0;
   int          nKnapCoverCuts   = 0;
   int          nMixIntRoundCuts = 0;
   int          nGomoryCuts      = 0;
   int          nTotalCuts       = 0;
   //---
   //--- this is typically coming from relaxed master problem
   //---  which has no defined integers (why not?) you are using
   //---  initalSolve, resolve which solves relaxation anyway...
   //---  can't the master and m_cutgenSI just stores the integers?
   //---
   //--- if you set integer on the fly and then unset
   //---   will that break things?
   //--- if master, need to set, if cutgenSI, don't...
   //---
   //int nInts = static_cast<int>(integerVars.size());
   //if(nInts > 0)
   // si->setInteger(&integerVars[0], nInts);
   //---
   //--- some CGLs need row activities too
   //---   currently, no easy way to set this
   //---
   OsiClpSolverInterface* cutGenClpSI =
      dynamic_cast<OsiClpSolverInterface*>(cutGenSI);
   assert(cutGenClpSI);
   const int nRows = cutGenClpSI->getNumRows();
   //---
   //--- calculate activity
   //---
   //TODO: note, this design never does cuts on cuts
   const CoinPackedMatrix* M   = cutGenClpSI->getMatrixByRow();
   double*                  act = new double[nRows];
   assert(act);
   M->times(xhat, act);
   //---
   //--- set primal column solution
   //---
   cutGenClpSI->setColSolution(xhat);
   //---
   //--- set primal row solution (i.e., activities)
   //---
   //write a si->setRowPrice for OsiClp
   //si->setRowPrice(act);//BAD NAME!
   bool            mustDeleteWS = true;
   CoinWarmStart* warmStart    = NULL;

   //TODO: check on crossover code - some speedups possible
   //  with a version that accepts memory - so not alloc/free
   //  too often

   switch (m_algo) {
   case CUT:
      //---
      //--- set master warm start in cgl SI
      //---
      warmStart = masterSI->getPointerToWarmStart(mustDeleteWS);
      cutGenClpSI->setWarmStart(warmStart);
      break;

   case PRICE_AND_CUT:
   case RELAX_AND_CUT:
      if (m_genGomory) {
         //---
         //--- crossover from xhat to basic solution
         //---
         //---
         //--- OsiClp::crossover
         //---   options - 0 no presolve (use primal and dual)
         //---             1 presolve (just use primal)
         //---             2 no presolve (just use primal)
         //---   basis   - 0 use all slack basis
         //---             1 try and put some in basis
         //---
#ifdef TRUNK_BUILD
         int crossOptions = 2;
         int crossBasis   = 1;
         //add obj cut
         //   obj >= dw-master obj - then generate gomory?
         //    or do that all the time for cuts?
         cutGenClpSI->crossover(crossOptions, crossBasis);
#endif
         //cutGenClpSI->resolve();//need?
         ///////////STOP -> getting all kinds of not violated
         ///// maybe try options=2, no presolve as it might be
         ///// screwing up the model?
         break;
      }

   default:
      break;
   }

   /*#ifdef __DECOMP_IP_CPX__
   OsiCpxSolverInterface * subprobSI_Cpx
      = dynamic_cast<OsiCpxSolverInterface*>(si);
   CPXENVptr cpxEnv = subprobSI_Cpx->getEnvironmentPtr();
   CPXLPptr  cpxLp  = subprobSI_Cpx->getLpPtr();

   int err = CPXcopystart( cpxEnv, cpxLp, NULL, NULL,
   		   const_cast<double*>( xhat ),
   		   const_cast<double*>( act  ),
   		   NULL, NULL );
   printf("Err=%d\n",err);fflush(stdout);
   assert(!err);
   #endif*/

   if (m_genClique) {
      UTIL_MSG(m_logLevel, 3,
               (*m_logStream) << "Calling cut generator: cliques\n";
              );
      m_genClique->generateCuts(*cutGenClpSI, osiCuts);
      nCliqueCuts = osiCuts.sizeCuts() - nTotalCuts;
      nTotalCuts  = osiCuts.sizeCuts();
   }

   if (m_genOddHole) {
      UTIL_MSG(m_logLevel, 3,
               (*m_logStream) << "Calling cut generator: cliques\n";
              );
      m_genOddHole->generateCuts(*cutGenClpSI, osiCuts);
      nOddHoleCuts = osiCuts.sizeCuts() - nTotalCuts;
      nTotalCuts   = osiCuts.sizeCuts();
   }

   if (m_genFlowCover) {
      UTIL_MSG(m_logLevel, 3,
               (*m_logStream) << "Calling cut generator: flow-covers\n";
              );
      m_genFlowCover->generateCuts(*cutGenClpSI, osiCuts);
      nFlowCoverCuts = osiCuts.sizeCuts() - nTotalCuts;
      nTotalCuts     = osiCuts.sizeCuts();
   }

   if (m_genKnapCover) {
      UTIL_MSG(m_logLevel, 3,
               (*m_logStream) << "Calling cut generator: knap-covers\n";
              );
      m_genKnapCover->generateCuts(*cutGenClpSI, osiCuts);
      nKnapCoverCuts = osiCuts.sizeCuts() - nTotalCuts;
      nTotalCuts     = osiCuts.sizeCuts();
   }

   if (m_genMixIntRound) {
      UTIL_MSG(m_logLevel, 3,
               (*m_logStream) << "Calling cut generator: mixint-round\n";
              );
      m_genMixIntRound->generateCuts(*cutGenClpSI, osiCuts);
      nMixIntRoundCuts = osiCuts.sizeCuts() - nTotalCuts;
      nTotalCuts       = osiCuts.sizeCuts();
   }

   if (m_genGomory) {
      UTIL_MSG(m_logLevel, 3,
               (*m_logStream) << "Calling cut generator: gomory\n";
              );
      m_genGomory->generateCuts(*cutGenClpSI, osiCuts);
      nGomoryCuts = osiCuts.sizeCuts() - nTotalCuts;
      nTotalCuts  = osiCuts.sizeCuts();
   }

   UTIL_MSG(m_logLevel, 3,
            (*m_logStream)
            << "Num clique     cuts= " << nCliqueCuts      << endl
            << "Num odd-hole   cuts= " << nOddHoleCuts     << endl
            << "Num flow-cover cuts= " << nFlowCoverCuts   << endl
            << "Num knap-cover cuts= " << nKnapCoverCuts   << endl
            << "Num mixed-int  cuts= " << nMixIntRoundCuts << endl
            << "Num gomory     cuts= " << nGomoryCuts      << endl;
           );
   //osiCuts.printCuts();
   int i;

   for (i = 0; i < osiCuts.sizeRowCuts(); i++) {
      CoinAssertDebug(osiCuts.rowCut(i).consistent());
      CoinAssertDebug(osiCuts.rowCut(i).consistent(*cutGenClpSI));
      CoinAssertDebug(!osiCuts.rowCut(i).infeasible(*cutGenClpSI));
      //CoinAssertDebug(osiCuts.rowCut(i).violated(xhat) > 1.e-5);
      DecompCutOsi* decompCut = new DecompCutOsi(osiCuts.rowCut(i));

      if (osiCuts.rowCut(i).violated(xhat) < DecompEpsilon) {
         UTIL_DEBUG(m_logLevel, 3,
                    (*m_logStream) <<
                    "WARNING: CGL cut " << i << " not violated." << endl;
                    osiCuts.rowCut(i).print();
                   );
      } else {
         newCuts.push_back(decompCut);
      }
   }

   UTIL_DEBUG(m_logLevel, 3,

   if (osiCuts.sizeColCuts() > 0) {
   (*m_logStream)
            << "WARNING: " << osiCuts.sizeColCuts()
            << " CGL col cuts found." << endl;
   }
             );

   //   if(nInts > 0)
   // si->setContinuous(&integerVars[0], nInts);
   if (mustDeleteWS && warmStart) {
      UTIL_DELPTR(warmStart);
   }

   UTIL_DELARR(act);
   return status;
}
