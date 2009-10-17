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

/*
  APP (no defaults):
  - createModel
  - solveSubproblem - not true if user defines partition (solve as IP)

  naming? solveSubproblem - but can be price or cut
  solveRelaxed in app?

  TODO: compressRow, compressCols
*/

/*
  really should be:

  C, PC branch to simplex / interior point
  RC    branch to subgradient / volume  

  OSI should NOT be in base class of DecompAlgo... 
  
  just wrap it everywhere... and only use OSI
  and the level where we call simplex or IPM (CLP has one)... try it!

  then, RC not needing OSI should work out just fine... 

  Get rid of all OSI in base!... Is that going to be possible?
  finish RC, then revert back and see how to handle... 

  get some stats on pi_i - pi_i+1 from iter to iter for stable
  vs unstable... how display this? norm
  
*/


#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
//#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "OsiNullSolverInterface.hpp"

#include "Decomp.h"
#include "DecompApp.h"
#include "DecompAlgo.h"
#include "DecompAlgoD.h"
#include "DecompAlgoC.h"

#include "DecompVarPool.h"
#include "DecompCutPool.h"

#include "DecompCutOsi.h"

#include "UtilHash.h"

#define PARM 0
//#define CUT_FIRST
//#define DO_INTERIOR


#ifdef __DECOMP_IP_CBC__
#include "CbcModel.hpp"
#include "CbcStrategy.hpp"
#endif


// --------------------------------------------------------------------- //
const char * DecompAlgo::m_classTag         = "\nD-ALGO   : ";

// --------------------------------------------------------------------- //
void DecompAlgo::setMasterBounds(const double * lbs,
                                 const double * ubs){
   int c;
   const int n_cols = m_masterSI->getNumCols();
   //keep in memory
   int    * index  = new int[n_cols];
   double * bounds = new double[2*n_cols];
   for(c = 0; c < n_cols; c++){
     index[c]      = c;
     bounds[2*c]   = lbs[c];
     bounds[2*c+1] = ubs[c];
   }
   //   for(c = 0; c < n_cols; c++){
   // //set all at once would be faster
   // m_masterSI->setColBounds(c, lbs[c], ubs[c]);
   //}
   m_masterSI->setColSetBounds(index, index+n_cols, bounds);
   UTIL_DELARR(index);
   UTIL_DELARR(bounds);
}

// --------------------------------------------------------------------- //
void DecompAlgo::tighten(int whichModel){
   //a chance to create an augmented model by running DecompAlgoC... 
   //this can be a derived Algo? or just an option for PC?
#ifdef CUT_FIRST
   //best bet is to create a new modelCore to work with
   //this means DecompAlgoC needs some way to transfer back a new model
   //almost like instead of a solve method, a "tighten" method
   DecompAlgoC CTighten(m_app);
   printf("\n\n\n========== CUT TIGHTEN START ==================\n\n\n");
   CTighten.m_isTightenAlgo = 1;
   CTighten.solve(whichModel);//think  
   m_whichCoreModel = whichModel;
   m_isTightenAlgo = 0;
   m_app->m_modelCore[whichModel]->nBaseRows = 
      m_app->m_modelCore[whichModel]->getNumRows();

   //???
   DecompVarList::iterator it;
   for(it = CTighten.m_vars.begin();
       it != CTighten.m_vars.end(); it++)
      m_initVars.push_back(*it);
   //appendVars
   //copy m_vars from CTigthen to 
   //m_initVars.insert(CTighten.m_vars.begin(),
   //                  CTighten.m_vars.end());
   printf("\nm_initVars.size() = %d", m_initVars.size());


   //actually - should send over z_CLB from CPM, so if we do this, it is 
   //only useful if the DW bound > z_CLB... and so we should use z_CLB as 
   //m_tlb? will it improve the bound? only if parts of polytope from sub
   //problems dominate the parts of polytope from cuts... 



   //this will tighten the basic model, but we want to then 
   //copy that into the model which will be actually used... 
   //copying seems silly... should be able to just point to it?
   printf("\n\n\n========== CUT TIGHTEN END ====================\n\n\n");
   //exit(1);
   //it is possible you will be IP feasible and done after this
   //need to return code back to main so user doesn't call solve after
   //tighten in that case... 

   //also - if any points are integral... they are candidates for columns
   //in PC... actually - ANY point at all is a candidate... right? ok - this
   //has now moved into the weird zone!

   //example - TSP, you start with 2-matching, and solve LP - many times
   //you get integral points, which are 2-matchings... then you add subtours...
   //but each of those x*'s coming from regular cutting plane steps are 
   //candndiate columns to start off DW - right? there is actually no need
   //to have the columns of DW be integral - of course, if not, then you 
   //kind of blow away the whole idea of DC... and integral separation
  
   //if you take all of those points to start off the DW columns, you are
   //actually guaranteed to get the cutting plane bound as DW master bound
   //which is what you want - rather than having to price out alot just to 
   //get there??

   int n_rows = m_app->m_modelCore[whichModel]->getNumRows();
   //this makes it so you have to run tighten before stabilized PC
   m_piEstimate = new double[n_rows];
   const double * pi = CTighten.m_masterSI->getRowPrice();
   memcpy(m_piEstimate, pi, n_rows * sizeof(double));
  
   for(int i = 0; i < n_rows; i++){
      printf("\nm_piEstimate[%d]: %g", i, m_piEstimate[i]);
   }

#endif
}


// --------------------------------------------------------------------- //
void DecompAlgo::initSetup(int whichModel){

   UTIL_MSG(m_app->m_param.LogLevel, 2,
	    (*m_osLog)
	    << m_classTag << "Initial Algo Setup"
	    << " (algo = " << DecompAlgoStr[m_algo] << ")";
	    );


   
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- initSetup()    ---- ";
              );

   //TODO: error code
   assert(whichModel <= static_cast<int>(m_app->m_modelCore.size()));
  
   //---
   //--- define the current application and model 
   //---
   m_whichModel = whichModel;  
   if(m_whichCoreModel == -1)
      m_whichCoreModel = whichModel;
#ifdef CUT_FIRST
   m_modelCore  = m_app->m_modelCore[m_whichCoreModel];
#else
   m_modelCore  = m_app->m_modelCore[m_whichModel];
#endif
   m_modelRelax = m_app->m_modelRelax[m_whichModel];
   assert(m_modelCore);
   assert(m_modelCore->M);

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
	      (*m_osLog)
	      << "\nModelCore  cols: " << m_modelCore->getNumCols()
	      << " rows: "             << m_modelCore->getNumRows();
	      if(m_modelRelax->M){
		 (*m_osLog)
		    << "\nModelRelax  cols: " << m_modelRelax->getNumCols()
		    << " rows: "             << m_modelRelax->getNumRows();
	      }
	      );
  
   //THINK:
   //for DC, we don't need the m_modelCore!, only m_modelRelax->.. and
   //if called from within DecompAlgoC, we dont' need that at all,
   //we can just use the m_m_modelRelax passed in?
  
  
   //STOP - specific to PC?? THINK??
   //.... but, also want for C, if doing DC!
   //.... so, do in any case, just check for empty?
   //---
   //--- create the solver interface for the relaxed problem [A',b] (c = 0)
   //---
   //THINK: this should only be done in the case we plan to augment
   //P' and/or if the user does not provide a solver for P'.... 
   //THINK: so this way, when we augment P' we have to augment osi
   //and m_m_modelRelax->M?? 
  
   //this is where we have to think - using OSI in base - is that correct?
  
   int n_integers;
   DecompConstraintSet * modelRelax = NULL;
  
   map<int, DecompConstraintSet*>::iterator mdit;  
   for(mdit = m_app->m_modelRelax.begin(); 
       mdit != m_app->m_modelRelax.end(); mdit++){
      modelRelax = mdit->second;    
      OsiSolverInterface * subprobSI = NULL;
      n_integers = static_cast<int>(modelRelax->integerVars.size());
      if(modelRelax->M){
         subprobSI = initSolverInterface();  
         subprobSI->loadProblem(*modelRelax->M, 
                                &modelRelax->colLB[0], 
                                &modelRelax->colUB[0], 0,
                                &modelRelax->rowLB[0],
                                &modelRelax->rowUB[0]);
         if(modelRelax->integerVars.size() > 0){
            subprobSI->setInteger(&modelRelax->integerVars[0],  n_integers);
         }
      }
      m_subprobSI.insert(make_pair(mdit->first, subprobSI));
   }
  

   //
   // --- open memory to store the current solution (in terms of x)
   //

   //think about this - how do you pass this info from algo to algo
   //if it is a member of this instantiation - ??
   //other algos you would not be calling solve() - you would just
   //be calling subfunctions - when get to hybrid, i think this will
   //make more sense
   m_xhat       = new double[m_modelCore->getNumCols()];
   //m_piEstimate = new double[m_modelCore->getNumRows()];

   //init from "tighten" .... or if not called, at least init from solving
   //original LP over A' and A''? -> just take duals of A''... hmmmm...
   //DecompAlgoC do over A' and A'' ?? think... 
  
   //UtilFillN(m_piEstimate, m_modelCore->getNumRows(), 0.0);

   //
   //--- PC: create an initial set of points F'[0] subseteq F' (c    + eps)
   //--- DC: create an initial set of points F'[0] subseteq F' (xhat + eps)
   //--- RC: do nothing - DecompAlgo base?? WHY - need an shat to get going
   //---  C: do nothing - DecompAlgo base
   //
   //THINK
   //if(m_algo != RELAX_AND_CUT)
   m_varsThisCall += generateInitVars(m_initVars);

   //THINK
#if 0
   if(m_algo == RELAX_AND_CUT){
      //set the reduced cost (u=0) for first iteration
      //(c - uA")s
      DecompVarList::iterator it;
      for(it = initVars.begin(); it != initVars.end(); it++){
         (*it)->setReducedCost((*it)->m_s.dotProduct(m_m_modelCore->objCoeff));
      }
   }
#endif
  
  
   //?? duMerle stablized version of PC master
   //
   //--- initialize solver interface (using OSI)
   //--- PC: min c(s lam)
   //---     A''s   lam   >= b, 
   //---     sum{s} lam_s  = 1,
   //---            lam_s >= 0, s in F'[0]
   //
   //TODO: this comment is in wrong place... specific to PC
   m_masterSI = initSolverInterface();
   createMasterProblem(m_initVars);

   m_auxMemPool.allocateMemory(m_modelCore->getNumCols(),
                               m_modelCore->getNumRows());

   //set timer limit
   m_stats.timerOverall.setLimit(m_app->m_param.LimitTime);
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- initSetup()    ----> ";
              );
}


// --------------------------------------------------------------------- //
DecompStat DecompAlgo::processNode(const int nodeIndex){  
   DecompPhase   phase              = PHASE_INIT;
   DecompStat    stat               = STAT_UNKNOWN;
   double        mostNegReducedCost = 0.0;

   UTIL_MSG(m_app->m_param.LogLevel, 2,
	    (*m_osLog)
	    << m_classTag << "Process Node " << nodeIndex
	    << " (algo = " << DecompAlgoStr[m_algo] 
	    << ", model = " << m_whichModel << ")";
	    ;
	    );

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- processNode() ---- ";
              );

   m_stats.timerDecomp.restart();
  
   //TODO: stats object - node stats, full solve stats
   m_nodeIndex       = nodeIndex;
   m_cutsThisRound   = 0;
   m_varsThisRound   = 0;
   m_cutsThisCall    = 0;
   m_varsThisCall    = 0;
   m_cutCallsTotal   = 0;
   m_priceCallsTotal = 0;
   m_cutCallsRound   = 0;
   m_priceCallsRound = 0;
   
   m_xhatIPBest      = 0;
   //m_tlb is the local (for this node) lb
   //m_tub is the global ub
   //we don't store the global lb.... 
   m_tlb             = -DecompInf; //? at least can set to global lb?
  
   // ---
   // --- initial update of the solution (dual and/or primal)
   // --- PC: take PARM steps of simplex
   // --- ?? DC: take PARM steps of simplex (INF case?)
   // --- RC: take PARM steps of subgradient
   // --- VC: take PARM steps of volume
   // ---
   if(m_algo != RELAX_AND_CUT)
      stat = solutionUpdate(phase, 99999, 99999);
   //setTrueLowerBound(0.0);
   if(m_algo == CUT){
      m_tlb = m_masterSI->getObjValue();
   }


   //don't want to do this if INF
   //STAT_INFEASIBLE means LP infeasible here... 
   //THINK: NO m_maseterSI in RC

   //STOP rewrite... 
   //printf("\nnumcols: %d", m_modelCore->getNumCols());
   if(stat != STAT_INFEASIBLE)
      recomposeSolution(m_masterSI->getColSolution(), m_xhat);
  
   //?? what about the INF case?? artificial columns? DC-ABCC version

   // ---
   // --- update the phase
   // ---
   if(m_algo != RELAX_AND_CUT)
      phase = phaseUpdate(phase, stat);
   else
      phase = PHASE_PRICE;
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 50,
              printCurrentProblem(m_masterSI, "masterProb", m_nodeIndex, 
                                  m_cutCallsTotal, m_priceCallsTotal);
              );
  
   while(phase != PHASE_DONE){
    
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << "\n\n========================= PHASE: " 
                 << DecompPhaseStr[phase] 
                 << "\tcutCalls: " << m_cutCallsRound
                 << " ( " << m_cutCallsTotal << " )"
                 << "\tpriceCalls: " << m_priceCallsRound
                 << " ( " << m_priceCallsTotal << " )"
                 << " masterObj: "
                 << m_masterSI->getObjValue()
                 << "\n";
                 );

      //---
      //--- check if we have exceeded time
      //---
      if(m_stats.timerOverall.isPast(m_app->m_param.LimitTime)){
         (*m_osLog) << "\nSTOP: Time Limit Hit\n";
         phase = PHASE_DONE;
      }
    
    
    
      DecompVarList newVars;
      DecompCutList newCuts;
    
      //these are from last time you called that func
      //these are updated in phase update
      //n_newVars = 0;
      //m_cutsThisCall = 0;
    
      //
      // --- calculate an upper bound 
      //


      //?? assumption, dual information is always available
      //
      // --- if dual information is available
      // --- PC: solve zSP(reduced cost) 
      // --- DC: solve zSP(reduced cost) ?? dual ray 
      // --- RC: do nothing
      //
      switch(phase){
      case PHASE_PRICE:
         //n_priceCalls++;
         //n_priceCallsTotal++;

         m_priceCallsRound++;
         m_priceCallsTotal++;
         //m_cutsThisRound = 0;
      
         // ---
         // --- after adding some rows, the columns in the var pool
         // --- might no longer be valid, so we need to re-expand everything
         // ---

         //this cannot work for DC
         if(!m_varpool.colsAreValid())
            m_varpool.reExpand(*m_modelCore, m_app->m_param.TolZero);

         //if RC, then vector u, no masterSI... 

         if(stat == STAT_FEASIBLE)
            m_varpool.setReducedCosts(m_masterSI->getRowPrice(), stat);
         else{
            //if doing RC, never called??
            const double * u = getDualRays(1)[0];
            m_varpool.setReducedCosts(u, stat);
            UTIL_DELARR(u);
         }
      
         // ---
         // --- attempt to generate some new variables with rc < 0
         // ---
         mostNegReducedCost = 0.0;
         //need to check for dups inside generateVars, since this 
         //flags whether or not to go for another round... if wait until 
         //pool - then that value is wrong, and phaseUpdate will be off?
         m_varsThisCall = generateVars(stat, newVars, mostNegReducedCost);
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                    (*m_osLog) << "\nm_varsThisCall: " << m_varsThisCall;
                    );
      
         m_varsThisRound += m_varsThisCall;

         //TODO: can stop DW as soon as lb is higher than
         //current ub.... like dual OBj limit

         //---
         //--- update the true lower bound
         //---
         //WRONG: LB = rc + alpha + ub'' (NEED ALPHA BACK IN)
         //this is wrong... AFTER A CUT, YOU HAVE TO PRICE BEFORE GET VALID BOUND
         
	 setTrueLowerBound(mostNegReducedCost); //what about for C?
         //NOTE: THIS CHECK IS NOT TRUE IN TREE!    
         //printf("\nm_tlb: %g, m_tub: %g, m_bestUpperBound: %g", 
         //       m_tlb, m_tub, m_bestUpperBound);

         if(m_nodeIndex == 0){
            assert(m_tlb <= (m_bestUpperBound + 1.0e-8));
            assert(m_tlb <= (m_tub            + 1.0e-8));
         }

         if(m_varsThisCall == 0){
            // ---
            // --- we have priced out, we are done with this phase
            // ---    
            //why do this here?
            //if(n_cutCallsTotal >= m_app->m_param.LimitTotalCutIters)
            //  phase = PHASE_DONE;
         }
         else{
        
            //THINK!!  should use var pool for volume etc ideas...
            if(m_algo != RELAX_AND_CUT){
               //--- 
               //--- add the newly generated variables to the var pool
               //---
               addVarsToPool(newVars);
          
               //---
               //--- add variables from the variable pool to the master problem
               //---
               addVarsFromPool();
            }
        
         }
         break;
      case PHASE_CUT:
         //n_cutCalls++;
         //n_cutCallsTotal++;

         m_cutCallsRound++;
         m_cutCallsTotal++;
         //m_varsThisRound = 0;
      

         // ---
         // --- after adding some cols, the rows in the cut pool
         // --- might no longer be valid, so we need to re-expand everything
         // ---
         //ok for PC, but if doing C, don't need right!?
         //in C, should be ok, should never be called as rowsAreValid = true
         if(!m_cutpool.rowsAreValid())
            m_cutpool.reExpand(m_vars, m_modelCore->getNumCols());
      
      
         //THINK: here is where you will do sep of xhat vs shat
         m_cutpool.calcViolations(m_xhat);

         //TODO: check if IP feasible - are we done? this part
         //might also find cuts... 
           
         // ---
         // --- attempt to generate some new cuts with vio > 0
         // ---
         m_cutsThisCall = generateCuts(newCuts);
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                    (*m_osLog) << "\nm_cutsThisCall: " << m_cutsThisCall;
                    );

         m_cutsThisRound += m_cutsThisCall;

         if(m_cutsThisCall == 0){
            //---
            //--- we have found no violated cuts, we are done with this phase
            //---
            //if(n_priceCallsTotal >= m_app->m_param.LimitTotalPriceIters)
            //  phase = PHASE_DONE;
         }
         else{
        
            //--- 
            //--- add the newly generated cuts to the cut pool
            //---
            //m_shat?
            addCutsToPool(m_xhat, newCuts, m_cutsThisCall);

            //---
            //--- add cuts from the cut pool to the master problem
            //---
            addCutsFromPool();        
         }
         break;
      case PHASE_DONE:
         break;
      default:
         assert(0);
      }

      if(phase != PHASE_DONE){
         // ---
         // --- perform a solution update
         // --- PC: take PARM steps of simplex
         // --- ?? DC: take PARM steps of simplex (INF case?)
         // --- RC: take PARM steps of subgradient
         // --- VC: take PARM steps of volume
         // ---
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                    (*m_osLog) << "\nm_cutsThisCall: " << m_cutsThisCall 
                    << " m_varsThisCall: " << m_varsThisCall;
                    );
         //THINK about this... 
         //price, cut (found cut), cut (no cuts found), why are we reoptimizing?
         if((m_cutsThisCall + m_varsThisCall) > 0){
            if(m_algo != RELAX_AND_CUT){
               UTIL_DEBUG(m_app->m_param.LogDebugLevel, 50,
                          printCurrentProblem(m_masterSI, "masterProb", m_nodeIndex, 
                                              m_cutCallsTotal, m_priceCallsTotal);
                          );
            }
            stat = solutionUpdate(phase, 99999, 99999);
         }

         //---
         //--- check if IP feasible (are we done?) 
         //--- TODO: for nonexplicity, also check user app isfeasible
         //---
         if(stat != STAT_INFEASIBLE){
            recomposeSolution(m_masterSI->getColSolution(), m_xhat);
            if(m_algo != RELAX_AND_CUT)
               UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
                          m_app->printOriginalSolution(m_modelCore->getNumCols(),
                                                       m_xhat, m_osLog);
                          );
         
            //this is NOT true if still pricing... THINK
            if(m_algo == CUT && (phase == PHASE_CUT || phase == PHASE_INIT))
               m_tlb = m_masterSI->getObjValue();

            //think about this - perhaps this should be in phase update, not here
#if 0
            if((m_varsThisCall == 0) && 
               (m_cutCallsTotal > 0) &&
               (m_cutsThisCall == 0) && 
               isIPFeasible(m_xhat)){          
               printf("\nIP FEASIBLE, we are done?");
               phase = PHASE_DONE;
            }
#endif
            //TODO
         }
      
      
         // adjust var effectivness, compress cols?
      
         //
         // --- if primal information is available
         //
         // THINK: cmp to m_tlb? if we meet it, might as well
         // stop processing this node and return
         //can the heur return all? or should it only return the best
         //i think alps manages this ok... 

         if(stat != STAT_INFEASIBLE){

            //this is checked again in phase update...
            //first, check to see if LP solution is already ip and user feas
            bool ipFeasible;
            bool appFeasible;
            if(isIPFeasible(m_xhat)){               
               //printf("\nxhat is IP FEASIBLE");
               if(m_app->APPisUserFeasible(m_xhat,
                                           m_modelCore->getNumCols(),
                                           m_app->m_param.TolZero)){
                  //printf("\nxhat is APP FEASIBLE");
                  DecompSolution * decompSol
                     = new DecompSolution(m_modelCore->getNumCols(),
                                          m_xhat, m_masterSI->getObjValue());
                  m_xhatIPFeas.push_back(decompSol);
               }
            }
            
            heuristics(m_xhat, m_xhatIPFeas);
            //m_app->APPheuristics(m_xhat, m_xhatIPFeas);
            
            //write functor for this - probably want solution pool
            //double bestQuality = DecompInf;
            vector<DecompSolution*>::iterator vi;
            DecompSolution * viBest = NULL;
            for(vi = m_xhatIPFeas.begin(); vi != m_xhatIPFeas.end(); vi++){
               const DecompSolution * xhatIPFeas = *vi;
               //xhatIPFeas->print();
               if(isIPFeasible(xhatIPFeas->getValues())){
                  //printf("\nIP Feasible from Heur Quality = %g",
                  //       xhatIPFeas->getQuality());
                  if(xhatIPFeas->getQuality() < m_tub){
                     m_tub = xhatIPFeas->getQuality();
                     viBest = *vi;
                  }
               }         
            }
            if(viBest){
               //save the best
               m_xhatIPBest = viBest;
               //printf("\nThe best:");
               //m_xhatIPBest->print();
            }	
            
            //clear out the rest
            //UtilDeleteVectorPtr(m_xhatIPFeas);
         }


         // ---
         // --- update the phase
         // ---
         //printf("\nphase: %s, stat: %s, n_cutCalls: %d, n_priceCalls: %d",
         //       DecompPhaseStr[phase], DecompStatStr[stat], 
         //       m_cutCallsRound, m_priceCallsRound);
         
         if(phase != PHASE_DONE){
            phase = phaseUpdate(phase, stat); 
         }
                  
         //TODO - process them!
      }

      if(m_algo == CUT){
         m_tlb = m_masterSI->getObjValue();
      }
   } //while(phase != PHASE_DONE)


   //THINK - if we get a new IP feasible point during search,
   //we can stop search for LB early

   //need to check again, if we get ip feasible in first LP
   //but this will cause dups... if we also find above?
   if(m_xhatIPFeas.size() == 0 && stat != STAT_INFEASIBLE){
      //this is checked again in phase update...
      //first, check to see if LP solution is already ip and user feas
      bool ipFeasible;
      bool appFeasible;
      if(isIPFeasible(m_xhat)){               
         printf("\nxhat is IP FEASIBLE");
         if(m_app->APPisUserFeasible(m_xhat,
                                     m_modelCore->getNumCols(),
                                     m_app->m_param.TolZero)){
            printf("\nxhat is APP FEASIBLE");
            DecompSolution * decompSol
               = new DecompSolution(m_modelCore->getNumCols(),
                                    m_xhat, m_masterSI->getObjValue());
            m_xhatIPFeas.push_back(decompSol);
            m_xhatIPBest = decompSol;
         }
      }
   }
            

   

   //check to see if we got lucky, is m_xhat IP feasible -
   //why do this in BcpsDecompTreeNode?


   //if we have an ip feasible point, 

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- processNode() ----> ";
              );


   m_stats.thisDecomp.push_back(m_stats.timerDecomp.timeElapsed());

   return stat;
}


// ------------------------------------------------------------------------ //
int DecompAlgo::heuristics(const double            * xhat,
                           vector<DecompSolution*> & xhatIPFeas){
   
   //dumb idea of solving the current master as an IP... if feasible
   //that is an incumbent... too hard? put time limit

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- heuristics ---- ";
              );
   
   m_app->APPheuristics(m_xhat, m_xhatIPFeas);
			

#if 0
   //this makes no sense in CUT
   if(m_algo == PRICE_AND_CUT){
   //this is silly - just check each new column as it comes in
   //if you get lucky... 

   OsiIpSolverInterface * si = new OsiIpSolverInterface();
   
#ifdef __DECOMP_LP_CLP__
#ifdef __DECOMP_IP_CBC__
   si->getModelPtr()->messageHandler()->
      setLogLevel(m_app->m_param.LogLpLevel);
   si->getRealSolverPtr()->messageHandler()->
      setLogLevel(m_app->m_param.LogLpLevel);
#endif
#endif
 
   printf("\nSET MSG LEVEL = %d", m_app->m_param.LogLpLevel);
   si->messageHandler()->setLogLevel(m_app->m_param.LogLpLevel);
   
   //TODO: faster if byCol?
   si->loadProblem(*m_masterSI->getMatrixByRow(),
                   m_masterSI->getColLower(),
                   m_masterSI->getColUpper(),
                   m_masterSI->getObjCoefficients(),
                   m_masterSI->getRowLower(),
                   m_masterSI->getRowUpper());
   //TODO: don't push in trivial col-bounds?

   int i;   
   for(i = 0; i < si->getNumCols(); i++){
      si->setInteger(i);
   }
   
   //si->writeMps("heurIP");
   //si->writeLp("heurIP");
   si->branchAndBound();

   const double * solution = si->getColSolution();
   for(i = 0; i < si->getNumCols(); i++){
      if(!UtilIsZero(solution[i]))
         printf("\nx[ %d ]: %g", i, solution[i]);
   }

   //does this mean integer optimal?
   if(si->isProvenOptimal()){
      //ok to use m_xhat? NO NO
      //recomposeSolution(solution, m_xhat);
      double * xhat = new double[m_modelCore->getNumCols()];
      recomposeSolution(solution, xhat);
      xhatIPFeas.push_back(new DecompSolution(m_modelCore->getNumCols(),
                                              xhat,
                                              si->getObjValue()));
      UTIL_DELARR(xhat);
   }
   UTIL_DELPTR(si); 
   //exit(1);
   }
#endif
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- heuristics()     ---->";
              );
   return 1;
}

// ------------------------------------------------------------------------ //
//this seems ok for C and PC... but what when we want to do DC within C THINK
int DecompAlgo::generateCuts(DecompCutList & newCuts){
   //pass in m_xhat??

   //TODO: cut shat
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- generateCuts() ---- ";
              );

   m_stats.timerOther1.restart();

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
              for(int i = 0; i < m_modelCore->getNumCols(); i++){
                 if(!UtilIsZero(m_xhat[i]))
                    (*m_osLog) << "\nxhat[" << i << "] : " << m_xhat[i];
              }
              );
             
   //here is where you ask the user how to generate cuts on x
   //shat??

#if 1
   if(m_app->m_param.CutDC){
      printf("\n\n==================================== IN DECOMP");
      DecompAlgoD D(m_app, m_xhat, m_modelCore->getNumCols());
    
      //either returns a set of cuts or a decomposition, have that wrap solve()?
      D.solveD(&newCuts);
      //who deletes this memory? better to pass in newCuts.. 
      printf("\n\n====================================OUT DECOMP");
   }
#endif
   //exit(1);

   m_app->generateCuts(m_xhat,
                       *m_modelCore,
                       *m_modelRelax,
                       newCuts);


   //---
   //--- attempt to generate CGL cuts on x??
   //--- the only way this is going to work, is if you carry
   //--- around another OSI problem instance completely in terms of x
   //--- only allow CGL to generate cuts on the original full formulation??
   //--- otherwise, we'd be allowing cuts on cuts - would have to add cuts
   //--- to master and to this version... hmmm... ughh

   //--- for some problems, you can't have the full original formulation
   //--- for PC you probably don't want it anyway... P' comes in from ep's
   //--- you can just generate cuts on Q"? that cut off xhat

   //--- but for C, you need Q' and Q"

   //--- m_masterSI holds the problem in terms of lambda (over Q")
   //--- m_subprobSI holds the problem in terms of x (over P')
   //---


   //TODO: need a switch user can use CGL or not... 

   //cannot do gomory cuts without tableau
   //THINK HERE!!

   //for C, this step is definitely NOT neccessary, 
   //just cut on masterSI... plus integers
  
   if(m_app->m_param.CutCGL){
      int i,j;

      //clearly, we don't want to do this for CPM
      //we want to use the actual Osi
      //this is expensive and wasteful?

      OsiSolverInterface * siCgl = new OsiNullSolverInterface();
    
      CoinAssertDebug(m_modelCore->M->getNumRows() 
                      == static_cast<int>(m_modelCore->rowLB.size()));
      CoinAssertDebug(m_modelCore->M->getNumRows() 
                      == static_cast<int>(m_modelCore->rowUB.size()));
      if(m_modelRelax->M){
         CoinAssertDebug(m_modelRelax->M->getNumRows() 
                         == static_cast<int>(m_modelRelax->rowLB.size()));
         CoinAssertDebug(m_modelRelax->M->getNumRows() 
                         == static_cast<int>(m_modelRelax->rowUB.size()));
      }
    
      char * colType = new char[m_modelRelax->getNumCols()];
      UtilFillN(colType, m_modelRelax->getNumCols(), 'C');
      for(unsigned int ii = 0; ii < m_modelRelax->integerVars.size(); ii++){
         //TODO: B vs I?
         colType[m_modelRelax->integerVars[ii]] = 'B';
      }


      //Q" plus cuts so far
      //ugh - think
      CoinPackedMatrix * Mcol = new CoinPackedMatrix;
      Mcol->reverseOrderedCopyOf(*m_modelCore->M);

      dynamic_cast<OsiNullSolverInterface*>(siCgl)
         ->loadDataAndSolution(
                               *m_modelCore->M,  
                               *Mcol,
                               &m_modelCore->colLB[0], 
                               &m_modelCore->colUB[0],                      
                               m_app->m_model.objCoeff, //think
                               &m_modelCore->rowLB[0], 
                               &m_modelCore->rowUB[0],
                               colType, m_xhat, DBL_MAX                          
                               );
    
      //or, could use bottomAppendPackedMatrix
      //THINK how to do
#if 0 
      if(m_modelRelax->M){
         for(i = 0; i < m_modelRelax->getNumRows(); i++){
            siCgl->addRow(m_modelRelax->M->getVector(i),
                          m_modelRelax->rowLB[i],
                          m_modelRelax->rowUB[i]);               
         }
         //also add in Q'?? also, cannot do this here everytime!!
         //save it.. THINK
      }

#endif

      //in case m_modelRelax is empty... this is wrong
      //actually, might want to just have integers as part of
      //overall model, not as one of the constraint sets?
      //siCgl->setInteger(&m_modelRelax->integerVars[0], 
      //                 static_cast<int>(m_modelRelax->integerVars.size()));

      //set col type instead - have model object hold colType?
      //dynamic_cast<OsiNullSolverInterface*>(siCgl)->setColType(colType);    
      //siCgl->setColSolution(m_xhat);

      //THINK: for knap covers, need RC, why can't we get the dual vector
      //from the x-space and calculate RC, given the optimal dual solution
      //to the DW master? that is, why can't we just recombine to get u just 
      //like x
      //cut pool management - not too many, passes, etc?

      OsiCuts cs;

      //TODO: 
      //idea of somehow using CPX cuts here... like running
      //just root node? 
      //CglGomory                cglGomory;

      if(m_app->m_param.CutCglKnapC){
         CglKnapsackCover cglKnap;
         cglKnap.generateCuts(*siCgl, cs);
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << "\nAfter Knaps cs.sizeCuts()      = "
		    << cs.sizeCuts();
		    );
      }
      if(m_app->m_param.CutCglFlowC){
         CglFlowCover cglFlowCover;
         cglFlowCover.generateCuts(*siCgl, cs);
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << "\nAfter FlowCovers cs.sizeCuts() = "
		    << cs.sizeCuts();
		    );
      }
      if(m_app->m_param.CutCglClique){
         CglClique cglClique;
         cglClique.setStarCliqueReport(false);
         cglClique.setRowCliqueReport(false);
         cglClique.generateCuts(*siCgl, cs);
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << "\nAfter Cliques cs.sizeCuts()    = "
		    << cs.sizeCuts();
		    );
      }
      if(m_app->m_param.CutCglMir){
         CglMixedIntegerRounding2 cglMirs;
         cglMirs.generateCuts(*siCgl, cs);
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << "\nAfter MIRS cs.sizeCuts()       = "
		    << cs.sizeCuts();
		    );
      }
      if(0){
         CglGomory cglGomory;
         cglGomory.generateCuts(*siCgl, cs);
         (*m_osLog) << "\nAfter Gomory cs.sizeCuts()     = " << cs.sizeCuts();
      }

      //probe uses reduced costs also (we can set it?)
      //CglProbing cglProbing; //ok for PC?
      //why are non-violated cuts being produced in CGL????
  







    
#ifdef __DECOMP_LP_CLP__
      //TODO: NOT GOING TO WORK FOR PC
      //something wrong with Gomory? 10/25/06
      //cglGomory.generateCuts(*siCgl, cs);
      //(*m_osLog) << "\nAfter Gomory cs.sizeCuts()     = " << cs.sizeCuts();
      //cglProbing.generateCuts(*siCgl, cs);
      //(*m_osLog) << "\nAfter Probing cs.sizeCuts()    = " << cs.sizeCuts();
#endif
      //exit(1);
    
#if 0
      //---
      //--- since we know xhat is LP feasible, crossover is simple
      //--- BUT... xhat is only LP feasible to Q"... NOT neccessarily Q'
      //--- in the context of DW...?? is that true? think... its a convex
      //--- combination of membesr of F'... so it will satisfy Q'.. well, 
      //--- not in full actually... depending on how many columns we have 
      //--- brough in so far... it is actually feasible for conv(E') inter Q"
      //---
      //--- so, if you were only doing over Q", this might work... 
      CoinWarmStartBasis ws;
      ws.setSize(m_modelCore->getNumCols(), 
                 m_modelCore->getNumRows() + m_modelRelax->getNumRows());
    
      double activity;
      int n_basic = 0;
      int m       = m_modelCore->getNumRows() + m_modelRelax->getNumRows();
      for(i = 0; i < m_modelCore->getNumCols(); i++){
         if(UtilIsZero(m_xhat[i] - m_modelCore->colLB[i])){
            ws.setStructStatus(i, CoinWarmStartBasis::atLowerBound);
         } else if(UtilIsZero(m_xhat[i] - m_modelCore->colUB[i])){
            ws.setStructStatus(i, CoinWarmStartBasis::atUpperBound);
         }
         else{
            ws.setStructStatus(i, CoinWarmStartBasis::basic);
            n_basic++;
         }
      }

      for(i = 0; i < m_modelCore->getNumRows(); i++){
         activity = m_modelCore->M->getVector(i).dotProduct(m_xhat); 
         printf("\ncore row: %d, act: %g, lb: %g, ub: %g, atLB? %d, atUB? %d",
                i, activity, m_modelCore->rowLB[i], m_modelCore->rowUB[i],
                UtilIsZero(activity - m_modelCore->rowLB[i]),
                UtilIsZero(activity - m_modelCore->rowUB[i]));
         if(UtilIsZero(activity - m_modelCore->rowLB[i])){
            ws.setArtifStatus(i, CoinWarmStartBasis::atLowerBound);
         } else if(UtilIsZero(activity - m_modelCore->rowUB[i])){
            ws.setArtifStatus(i, CoinWarmStartBasis::atUpperBound);
         }
         else{
            printf(" -> basic");
            ws.setArtifStatus(i, CoinWarmStartBasis::basic);
            n_basic++;
         }
      }
      for(i = 0; i < m_modelRelax->getNumRows(); i++){
         activity = m_modelRelax->M->getVector(i).dotProduct(m_xhat);  
         printf("\nrelax row: %d, act: %g, lb: %g, ub: %g",
                i, activity, m_modelRelax->rowLB[i], m_modelRelax->rowUB[i]);
         if(UtilIsZero(activity - m_modelRelax->rowLB[i])){
            ws.setArtifStatus(i + m_modelCore->getNumRows(), 
                              CoinWarmStartBasis::atLowerBound);
         } else if(UtilIsZero(activity - m_modelRelax->rowUB[i])){
            ws.setArtifStatus(i + m_modelCore->getNumRows(),  
                              CoinWarmStartBasis::atUpperBound);
         }
         else{
            printf(" -> basic");
            ws.setArtifStatus(i + m_modelCore->getNumRows(), 
                              CoinWarmStartBasis::basic);
            n_basic++;
         }
      }
      printf("\nm = %d, n_basic = %d", m, n_basic);
      ws.print();
      assert(m == n_basic);
    
    
      siCgl->setWarmStart(&ws);
    
      cglGomory.generateCuts(*siCgl, cs);
      (*m_osLog) << "\nAfter GOMORY cs.sizeCuts()     = " << cs.sizeCuts();
    
    
#endif
    
#if 0
      OsiSolverInterface * siCgl2 = new OsiXxxSolverInterface();
      siCgl2->loadProblem(*m_modelCore->M,  //Q" plus cuts so far
                          &m_modelCore->colLB[0], 
                          &m_modelCore->colUB[0], 
                          //&m_modelCore->objCoeff[0], 
                          m_modelCore->objCoeff, 
                          &m_modelCore->rowLB[0], 
                          &m_modelCore->rowUB[0]);
      siCgl2->setInteger(&m_modelRelax->integerVars[0], 
                         m_modelRelax->integerVars.size());
      siCgl2->setColSolution(m_xhat);
    
      CoinWarmStartBasis ws;
      ws.setSize(m_modelCore->getNumCols(), 
                 m_modelCore->getNumRows());
    
      double activity;
      int n_basic = 0;
      int m       = m_modelCore->getNumRows();
      for(i = 0; i < m_modelCore->getNumCols(); i++){
         if(UtilIsZero(m_xhat[i] - m_modelCore->colLB[i])){
            ws.setStructStatus(i, CoinWarmStartBasis::atLowerBound);
         } else if(UtilIsZero(m_xhat[i] - m_modelCore->colUB[i])){
            ws.setStructStatus(i, CoinWarmStartBasis::atUpperBound);
         }
         else{
            ws.setStructStatus(i, CoinWarmStartBasis::basic);
            n_basic++;
         }
      }
    
      for(i = 0; i < m_modelCore->getNumRows(); i++){
         activity = m_modelCore->M->getVector(i).dotProduct(m_xhat); 
         printf("\ncore row: %d, act: %g, lb: %g, ub: %g, atLB? %d, atUB? %d",
                i, activity, m_modelCore->rowLB[i], m_modelCore->rowUB[i],
                UtilIsZero(activity - m_modelCore->rowLB[i]),
                UtilIsZero(activity - m_modelCore->rowUB[i]));
         if(UtilIsZero(activity - m_modelCore->rowLB[i])){
            ws.setArtifStatus(i, CoinWarmStartBasis::atLowerBound);
         } else if(UtilIsZero(activity - m_modelCore->rowUB[i])){
            ws.setArtifStatus(i, CoinWarmStartBasis::atUpperBound);
         }
         else{
            printf(" -> basic");
            ws.setArtifStatus(i, CoinWarmStartBasis::basic);
            n_basic++;
         }
      }
      printf("\nm = %d, n_basic = %d", m, n_basic);
      ws.print();
      assert(m == n_basic);
    
      siCgl2->setWarmStart(&ws);
    
      cglGomory.generateCuts(*siCgl2, cs);
      (*m_osLog) << "\nAfter GOMORY cs.sizeCuts()     = " << cs.sizeCuts();
      UTIL_DELPTR(siCgl2);
#endif
    
    
      //why is CGL probing returning column cuts that are not violated?
      //even though, should we be chaning column bounds?
    
      //printf("\ncs cuts     : %d", cs.sizeCuts());
      //printf("\ncs row cuts : %d", cs.sizeRowCuts());
      for(i = 0; i < cs.sizeRowCuts(); i++){
         //TODO: good for debugging - create this as member of DECOMP cut
         //classes too?
         //TODO: need to distinguish between col cuts and row cuts?? in decomp
         //part?
         //need a cutpool for col cuts??
         CoinAssertDebug(cs.rowCut(i).consistent());
         CoinAssertDebug(cs.rowCut(i).consistent(*siCgl));
         CoinAssertDebug(!cs.rowCut(i).infeasible(*siCgl));
      
         DecompCutOsi * decompCut = new DecompCutOsi(cs.rowCut(i));
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
		    (*m_osLog) << "\nCS ROW CUT " << i;
		    decompCut->print(m_osLog);
		    );
         newCuts.push_back(decompCut);
      }
    
      //so that it works in DW land, we have to treat these
      //column cuts as regular cuts... 
      for(i = 0; i < cs.sizeColCuts(); i++){
         CoinAssertDebug(cs.colCut(i).consistent());
         CoinAssertDebug(cs.colCut(i).consistent(*siCgl));
         CoinAssertDebug(!cs.colCut(i).infeasible(*siCgl));
      
         OsiColCut        & cCut = cs.colCut(i);
         const CoinPackedVector & lbs  = cCut.lbs();
         const CoinPackedVector & ubs  = cCut.ubs();
      
         //x >= l
         if(lbs.getNumElements() > 0){
            const int    * ind = lbs.getIndices();
            const double * els = lbs.getElements();
            for(j = 0; j < lbs.getNumElements(); j++){      
               CoinPackedVector v;
               v.insert(ind[j], 1.0);
               OsiRowCut rc;
               rc.setRow(v);
               rc.setLb(els[j]);
               rc.setUb(DecompInf);
        
               DecompCutOsi * decompCut = new DecompCutOsi(rc);
	       UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
			  (*m_osLog) << "\nCS COL >= CUT " << j;
			  decompCut->print(m_osLog);
			  );
               newCuts.push_back(decompCut);
            }
         }
      
         //x <= u
         if(ubs.getNumElements() > 0){
            const int    * ind = ubs.getIndices();
            const double * els = ubs.getElements();
            for(j = 0; j < ubs.getNumElements(); j++){      
               CoinPackedVector v;
               v.insert(ind[j], 1.0);
               OsiRowCut rc;
               rc.setRow(v);
               rc.setLb(-DecompInf);
               rc.setUb(els[j]);
          
               DecompCutOsi * decompCut = new DecompCutOsi(rc);
	       UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
			  (*m_osLog) << "\nCS COL <= CUT " << j;
			  decompCut->print(m_osLog);
			  );
               newCuts.push_back(decompCut);
            }
         }
      }
    
      UTIL_DELPTR(siCgl);
      delete [] colType;//THINK
      delete Mcol;
   }
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- generateCuts() ----> ";
              );
  
   m_stats.thisGenCuts.push_back(m_stats.timerOther1.timeElapsed()); 


   //think

  
   return static_cast<int>(newCuts.size());
}


// ------------------------------------------------------------------------ //
int DecompAlgo::generateVars(const DecompStat   stat,
                             DecompVarList    & newVars, 
			     double           & mostNegReducedCost){
  
   // ---
   // --- solve min{s in F' | RC[s]} to generate new variables
   // --- 
   // --- if LP was feasible,   then RC[s] = c.s - (uhat A''s + alpha)
   // --- if LP was infeasible, then RC[s] = (uhat A''s + alpha), dual ray
   // ---
   // --- The master LP was formed in the following order
   // ---   (A''s) lam[s]  >= b'' - from the original core [A'', b'']
   // ---   sum{s} lam[s]   = 1   - convexity constraint
   // --- But, we may have added cuts - which conceptually are added
   // --- into [A'', b'']. But, in reality they are simply appended to
   // --- the end of the LP matrix. So, when we get back the dual vector, 
   // --- we have to be aware of where alpha actually is.
   // ---
   //

   //THINK: do we really have to do it this way? seems a bit wasteful
   //int whichModel;
   const int      m          = m_masterSI->getNumRows();
   const int      n_corecols = m_modelCore->getNumCols();
   const double * u          = NULL;
   double       * redCostX   = NULL;

   double alpha;

   //FEAS RC = c - (uA''s + alpha) -> (c-uA'')s - alpha
   //INF  RC =     (us    + alpha) -> (u)s      + alpha   : WEIRD

   m_stats.timerOther1.restart();

   //THINK
   int m_nBaseCoreRows = m_modelCore->nBaseRows;
   if(m_algo == DECOMP)
      m_nBaseCoreRows = m_masterSI->getNumRows() - 1;
  
   //THINK - different for RC!  
   if(stat == STAT_FEASIBLE){
      u     = m_masterSI->getRowPrice();
      alpha = -u[m_nBaseCoreRows];    //dual of convexity constraint
   }
   else{  
      //---
      //--- user must delete memory (dual ray) in this case
      //---
      u     = getDualRays(1)[0];   //OK - if PRESOLVE finds inf?? 
      //alpha = u[m_nBaseCoreRows];  //dual of convexity constraint
      alpha = u[m_nBaseCoreRows];  //dual of convexity constraint
   }
   //THINK, is dual ray = (uhat A'' s + alpha) already?

  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- generateVars() ---- ";
              );
  
 
   double * u_adjusted = new double[m-1];       //core + cuts (no convexity)
   //TODO - util function? or should we use Coin Helper funcs?
   CoinDisjointCopyN(u, m_nBaseCoreRows, u_adjusted);
  
   //printf("\nm_cuts.size              = %d", m_cuts.size());
   //printf("\nm_masterSI->getNumRows() = %d", m_masterSI->getNumRows());
   //printf("\nm_masterSI->getNumCols() = %d", m_masterSI->getNumCols());
   //printf("\nm_nBaseCoreRows          = %d",  m_nBaseCoreRows);
#ifndef EXPLICIT_BOUNDS
   assert(static_cast<int>(m_cuts.size()) == m - m_nBaseCoreRows - 1);
#endif
   CoinDisjointCopyN(u + m_nBaseCoreRows + 1, 
                     m - m_nBaseCoreRows - 1,
                     u_adjusted + m_nBaseCoreRows);
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
              for(int i = 0; i < m_nBaseCoreRows; i++){  
                 printf("\nu[%d]: %g", i, u[i]);
              }  
              printf("\nalpha = %g", alpha);
              for(int i = 0; i < m-1; i++){  
                 printf("\nu_adj[%d]: %g", i, u_adjusted[i]);
              }
              );

   //THINK
   //the reason we use m_modelCore->M versus masterSI, is because
   //in RC, masterSI is not updated? otherwise, what is the whole
   //point of carrying around m_modelCore??
  
   DecompVarList potentialVars;
   redCostX = new double[n_corecols]; // (c - uhatA") in x-space
   if(m_algo == DECOMP){
      memcpy(redCostX, u_adjusted, n_corecols * sizeof(double));

#if 0
      //THINK? is this right? and WHY!!
      m_masterSI->getMatrixByCol()->dumpMatrix();
      m_masterSI->getMatrixByRow()->dumpMatrix();
      const CoinPackedMatrix * Mrow = m_masterSI->getMatrixByRow();
#if 1
      const int * index_ = Mrow->getIndices();
      const double * element_ = Mrow->getElements();
      memset(redCostX, 0.0, n_corecols * sizeof(double));
      for (int i = n_corecols - 1; i >= 0; --i) {
         double y_i = 0;
         const int last = Mrow->getVectorLast(i);
         for (int j = Mrow->getVectorFirst(i); j < last; ++j){
            if(index_[j] == m_nBaseCoreRows)
               continue;
            y_i += u_adjusted[index_[j]] * element_[j];
         }
         redCostX[i] = y_i;
      }
#endif
#endif

      //problem - don't want convexity in this multiplication?
      //m_masterSI->getMatrixByRow()->transposeTimes(u, redCostX);
      //for(int i = 0; i < n_corecols; i++){    
      //  printf("\nDredCostX[%d]: %g", i, redCostX[i]);
      //}  
   }
   else{    
      //core = A'' (m x n), DW = A'' s
      //for DECOMP, m_masterSI->M = 1 s (matrix for DECOMP problem)
      m_modelCore->M->transposeTimes(u_adjusted, redCostX);
      //for(int i = 0; i < n_corecols; i++){      
      //  printf("\nobjCoeff[%d]: %g redCostX[%d]: %g",
      //         i, m_app->m_model.objCoeff[i], i, redCostX[i]);
      // }
      if(stat == STAT_FEASIBLE){
         for(int i = 0; i < n_corecols; i++){    
            redCostX[i] = m_app->m_model.objCoeff[i] - redCostX[i];
            // printf("\nobjCoeff[%d]: %g redCostX[%d]: %g",
            //       i, m_app->m_model.objCoeff[i], i, redCostX[i]);
         }
      }
      else{
         //TODO: multiple farkas
         //for(int i = 0; i < n_corecols; i++){
         //redCostX[i] = -redCostX[i];//??
         //}
      }
   }
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 50,
              printCurrentProblem(m_subprobSI[m_whichModel], "subProb",
                                  m_nodeIndex,
                                  m_cutCallsTotal,
                                  m_priceCallsTotal);
              );

   //all models?
   //TODO: stat return, restrict how many? pass that in to user?
   if(!m_app->m_param.PriceMultiPoly){
      solveRelaxed(m_whichModel,
                   redCostX, 
                   m_app->m_model.objCoeff, 
                   alpha, //really -alpha
                   n_corecols, 
                   true, 
                   true, 
                   m_subprobSI[m_whichModel], 
                   potentialVars);
   }
   else{    
      map<int, DecompConstraintSet*>::iterator mdit;    
      for(mdit = m_app->m_modelRelax.begin(); 
          mdit != m_app->m_modelRelax.end(); mdit++){
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 50,
                    (*m_osLog) << "Solving relaxation " << mdit->first;
                    );
         solveRelaxed(mdit->first,
                      redCostX, 
                      m_app->m_model.objCoeff, 
                      alpha, //really -alpha
                      n_corecols, 
                      true, 
                      true, 
                      m_subprobSI[mdit->first],
                      potentialVars);
      }
   }
    

   //TODO: deal with duplicate columns - use ideas, Yan/COIN, hash key
   //TODO: add to DecompVar some indicator label from which poly it came from

   //only take the best if doing multi; best=max?
   DecompVarList::iterator it;
   double varRedCost;


   /*double leastNegReducedCost = -999999;//assume one per poly
   //if(m_app->m_param.PriceMultiPoly){
      for(it = potentialVars.begin(); it != potentialVars.end(); it++){
	 varRedCost = (*it)->getReducedCost();
	 UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		    (*m_osLog) << m_classTag << "alpha:" << alpha 
		    << " varRedCost: " << varRedCost;
		    );
         if(varRedCost > leastNegReducedCost){
            leastNegReducedCost = varRedCost;      
	 }
      }
      mostNegReducedCost = leastNegReducedCost;
      if(mostNegReducedCost < -1.0e-8){	 
	 for(it = potentialVars.begin(); it != potentialVars.end(); it++){
	    varRedCost = (*it)->getReducedCost();
	    UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		       (*m_osLog) << m_classTag << "alpha:" << alpha 
		       << " varRedCost: " << varRedCost;
		       );
	    if(varRedCost >= mostNegReducedCost){
	       newVars.push_back(*it);
	    }
	    else{
	       UTIL_DELPTR(*it);
	    }
	 }
	 }*/
      //}
      //else{

   for(it = potentialVars.begin(); it != potentialVars.end(); it++){
      varRedCost = (*it)->getReducedCost();
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
		 (*m_osLog) << m_classTag << "alpha:" << alpha 
		 << " varRedCost: " << varRedCost;
		 );
      if(varRedCost < -1.0e-8){ //TODO: strict, -dualTOL?
	 newVars.push_back(*it);
         if(varRedCost < mostNegReducedCost)
            mostNegReducedCost = varRedCost;      
      }
      else{
	 UTIL_DELPTR(*it);
	 //THINK: have to delete this memory?
      }  
      // }
      
   }
   //TODO: CHANGE NAME, adding alpha so TrueLB is correct - BAD NAME!!
   //alpha is already added in during solveRelaxed?

   //z_ = rc + alpha + ub'' (this alpha is really -alpha!)
   //THINK
   //mostNegReducedCost -= alpha;
   potentialVars.clear(); //THINK? what does clear do exactly ?

   for(it = newVars.begin(); it != newVars.end(); it++)
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*it)->print(m_osLog, m_app);
                 );

   //---
   //--- free local memory
   //---
   UTIL_DELARR(u_adjusted);
   UTIL_DELARR(redCostX);
   if(stat != STAT_FEASIBLE){
      UTIL_DELARR(u);
   }


   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- generateVars() ----> ";
              );

   m_stats.thisGenVars.push_back(m_stats.timerOther1.timeElapsed()); 

   return static_cast<int>(newVars.size());
}

// ------------------------------------------------------------------------- //
OsiSolverInterface * DecompAlgo::initSolverInterface(){
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- initSolverInterface() ---- ";
              );
    
   OsiLpSolverInterface * si = new OsiLpSolverInterface();  
   //TODO: different for master and subproblem solver?
   //si->setDblParam(OsiDualTolerance,   m_app->m_param.dualTol);
   //si->setDblParam(OsiPrimalTolerance, m_app->m_param.primalTol);
   si->messageHandler()->setLogLevel(m_app->m_param.LogLpLevel);

   //THINK this is called for m_subprobSI and m_masterSI
   //so can't do this here, need two pointers? and do each time init OSI?
   //m_masterCLP = si->getModelPtr(); 

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- initSolverInterface() ---->";
              );
   return si;
}

// ------------------------------------------------------------------------- //
bool DecompAlgo::isIPFeasible(const double * x,
                              const double   feasTol,
                              const double   intTol){

   //what does this mean in context of RC?
   if(!isLPFeasible(x, feasTol)){
      //printf("\nERROR??? NOT LP FEASIBLE");
      return false;
   }
  
   //assume integers are only defined in m_modelRelax
   //what about the case the uesr wants the relaxd problem to not have integral
   //but wants to define integers for overall problem - should not be member of
   //relaxed? or both?
   int c;

   vector<int>::iterator vi;  
   for(vi = m_modelRelax->integerVars.begin(); 
       vi != m_modelRelax->integerVars.end(); vi++){    
      c = *vi;
      if(!UtilIsIntegral(x[c], intTol))
         return false;
   }
  

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              m_app->printOriginalSolution(m_modelCore->getNumCols(), x);
              )
  
      return true;
}

// ------------------------------------------------------------------------- //
bool DecompAlgo::isLPFeasible(const double * x,
                              const double   feasTol){
   //---
   //--- The base isFeasible assumes a full explicit description
   //--- in m_modelCore and m_modelRelax (plus integrality). The user
   //--- app can also define an isFeasible - which will be checked
   //--- first.. THINK: but, then can't really use bool
   //---
   //--- TODO: isFeasible for TSP, etc might also be able to return
   //--- cuts... make this possible
   //---

   //double * MCore_x  = new double[m_modelCore->getNumCols()];
   //double * MRelax_x = new double[m_modelCore->getNumCols()];

   //THINK: could use submatrixOf since only need base rows,
   //but, the implementation seems very slow... maybe better to 
   //just loop and do each row
   //m_modelCore->M->times(m_xhat, MCore_x);
   //m_modelCore->M->times(m_xhat, MRelax_x);
  
   //---
   //--- do we satisfy all column bounds
   //---
   int c;
   double obj = 0.0;
   const double * objCoeff = m_app->m_model.objCoeff;
   for(c = 0; c < m_modelCore->getNumCols(); c++){
      obj += x[c] * objCoeff[c];
      if(x[c] < (m_modelCore->colLB[c] - feasTol)){
         //(*m_osLog) << "\nCoreIsLpFeas x[" << c << "]: " << x[c];
         //if(m_modelCore->colLB[c] > -DecompInf)
         //   (*m_osLog) << " lb: " << m_modelCore->colLB[c];
         //(*m_osLog) << " lb: -INF";
         return false;
      }
      else if(x[c] > (m_modelCore->colUB[c] + feasTol)){
         //if(m_modelCore->colLB[c] > -DecompInf)
         //   (*m_osLog) << " ub: " << m_modelCore->colLB[c];
         //(*m_osLog) << " ub: INF";
         return false;      
      }
   }

   //---
   //--- do we satisfy all row bounds
   //---
   int r;
   double ax;
   //if M is coloredered, then this is wrong
   for(r = 0; r < m_modelCore->nBaseRows; r++){
      ax = m_modelCore->M->getVector(r).dotProduct(x);
      if(ax < (m_modelCore->rowLB[r] - feasTol)){
         //(*m_osLog) << "\nCoreIsLPFeas ax[" << r << "]: " << ax;
         //if(m_modelCore->rowLB[r] > -DecompInf)
         //   (*m_osLog) << " lb: " << m_modelCore->rowLB[r];
         //else
         //   (*m_osLog) << " lb: -INF";
         //(*m_osLog) << "\n";
         //UtilPrintPackedVector(m_modelCore->M->getVector(r), m_osLog, m_app);
         return false;
      }
      else if(ax > (m_modelCore->rowUB[r] + feasTol)){
         //(*m_osLog) << "\nCoreIsLPFeas ax[" << r << "]: " << ax;
         //if(m_modelCore->rowUB[r] < DecompInf)
         //   (*m_osLog) << " ub: " << m_modelCore->rowUB[r]; 
         //else
         //   (*m_osLog) << " ub: INF";
         //(*m_osLog) << "\n";
         //UtilPrintPackedVector(m_modelCore->M->getVector(r), m_osLog, m_app);
         return false;
      }
   }
   if(m_modelRelax->M){
      for(r = 0; r < m_modelRelax->getNumRows(); r++){
         ax = m_modelRelax->M->getVector(r).dotProduct(x);
         if(ax < (m_modelRelax->rowLB[r] - feasTol)){
            //(*m_osLog) << "\nRelaxIsLPFeas ax[" << r << "]: " << ax;
            //if(m_modelRelax->rowLB[r] > -DecompInf)
            //   (*m_osLog) << " lb: " << m_modelRelax->rowLB[r];
            //else
            //   (*m_osLog) << " lb: -INF";
            //(*m_osLog) << "\n";
            //UtilPrintPackedVector(m_modelRelax->M->getVector(r), 
	    //                     m_osLog, m_app);
            return false;
         }
         else if(ax > (m_modelRelax->rowUB[r] + feasTol)){
            //(*m_osLog) << "\nRelaxIsLPFeas ax[" << r << "]: " << ax;
            //if(m_modelRelax->rowUB[r] < DecompInf)
            //   (*m_osLog) << " ub: " << m_modelRelax->rowUB[r]; 
            //else
            //   (*m_osLog) << " ub: INF";
            //(*m_osLog) << "\n";
            //UtilPrintPackedVector(m_modelRelax->M->getVector(r), 
	    //                     m_osLog, m_app);
            return false;         
         }
      }
   }  
   return true;
}

// ------------------------------------------------------------------------- //
//either RC gets a new setTrueLowerBound, or have a base derivation
//of getRowPrice, where RC gets it differently... 
//what the heck does this do for C??
//z_ = rc + alpha + ub''
void DecompAlgo::setTrueLowerBound(const double   mostNegReducedCost){  




   //this is suppose to be:
   //   zhat + rc* <= z* <= zhat, where
   //   zhat is the LP obj of current DW master
   //   we have been calculating this as ub + alpha (because by
   //   duality, that should equal the current LP bound)... however,
   //   that duality argument will change if any columns have bounds
   //   which we have because of the way we are doing branching??
   






   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- setTrueLowerBound()   ---- ";
              );
   //do we want to do this in RC, no m_masterSI
   if(m_algo == RELAX_AND_CUT || m_masterSI->isProvenOptimal()){//THINK
      //THIS is WRONG again for PC... you'll have alpha in here if any cuts
      //got added...
      //double tlb = mostNegReducedCost + calcConstant(m_modelCore->getNumRows(),
      //						   getRowPrice());

      double tlb = mostNegReducedCost + m_masterSI->getObjValue();

      //(*m_osLog) << "\nmostNegReducedCost: " << mostNegReducedCost 
      //           << "\nBEFORE tlb: " << tlb << " m_tlb: " 
      //           << (m_tlb < -1.0e20 ? -1.0e20 : m_tlb);
    
      m_tlb = m_tlb < tlb ? tlb : m_tlb; //tlb will be monotonic (during PRICE)
    
      //(*m_osLog) << "\nAFTER tlb: " << tlb << " m_tlb: " 
      //           << (m_tlb < -1.0e20 ? -1.0e20 : m_tlb);
   }
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- setTrueLowerBound()   ---->";
              );
}

// ------------------------------------------------------------------------- //
double DecompAlgo::calcConstant(const int      m, 
                                const double * u){//this is ub
   //TODO: this should be a denseDotProduct util func
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- calcConstant()        ---- ";
              );
   //const double * rhs = getRightHandSide();

   //core holds A'' + cuts
   //but u is on A'' + convex + cuts

   int convexity_index = m_modelCore->nBaseRows;
   const double * rhs = &m_modelCore->rowRhs[0];
   double constant = 0.0;
   for(int i = 0; i < convexity_index; i++){
      constant += u[i] * rhs[i];    
      //(*m_osLog) << "\ni: " << i << " u: " << u[i] << " rhs: "
      //           << rhs[i] << " const: " << constant;    
   }
   for(int i = (convexity_index+1); i <= m; i++){
      constant += u[i] * rhs[i-1];    
      //(*m_osLog) << "\ni: " << i << " u: " << u[i] << " rhs: "
      //           << rhs[i-1] << " const: " << constant;    
   }
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- calcConstant()        ---->";
              );
   return constant;
}

// ------------------------------------------------------------------------- //
//member of varpool versus algo class? different for DC??
void DecompAlgo::addVarsToPool(DecompVarList & newVars){

   //TODO: how to get rid of duplicate vars - use some kind of hash map!!
   //what we have now is really expensive! will see in profile

   double * denseCol = NULL;
   CoinPackedVector * sparseCol = NULL;

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- addVarsToPool()       ---- ";
              );
  
   if(m_algo != DECOMP)
      denseCol = new double[m_modelCore->getNumRows() + 1];
  
   DecompVarList::iterator li;//const?
   for(li = newVars.begin(); li != newVars.end(); li++){
      //(*li)->print();
    
      // --- 
      // --- get dense column = A''s, append convexity constraint on end 
      // --- THINK: PC specific
      // ---
      //M contains A'', OSI does not - that contains lambda space
      //but in decomp_and_cut, we don't need m_modelCore->M
    
      //the column that eventually gets added to masterSI, masterSI
      //has A'' (original) + convexity + cuts, in that order
      //so the column must match that order
      //but m_modelCore->M is only A'' + cuts, no convexity constraint

      //so need something different for RC vs PC vs ... ?
    
      if(m_algo == DECOMP){//UGH - derive own!
         sparseCol = new CoinPackedVector((*li)->m_s);
         UtilPrintPackedVector((*li)->m_s);
         sparseCol->insert(m_numOrigCols, 1.0);//convexity 
         UtilPrintPackedVector(*sparseCol);
      }
      //    else 
      //if(m_algo == PRICE_AND_CUT){
      //if this is needed, then just derive it!
      //that is the point of having virtual decompalgo methods
      //printf("\nHERE, calc sparseCol");
      //const CoinPackedMatrix * App = m_masterSI->getMatrixByRow();
      //App->times((*li)->m_s, denseCol);

      // ---
      // --- creat a sparse column from the dense column
      // ---
      //sparseCol 
      //  = UtilPackedVectorFromDense(m_masterSI->getNumRows(),
      //                              denseCol, m_app->m_param.TolZero);
      //}
      else{
         //modelCore->M = A'', plus possible cuts, but without convexity
         //this is WRONG as soon as a cut is added!

         //this is slow -- use BLAS?

         //BLAS won't help

         //this is slow to use sparse m_s because timesMinor needs to do
         //an expensive lookup x[index[j]]
         //which is faster, make dense, then this or ...

         (*li)->fillDenseArr(m_modelCore->getNumCols(),
                             m_auxMemPool.dblArrNCoreCols);
         //m_modelCore->M->times((*li)->m_s, denseCol);
         m_modelCore->M->times(m_auxMemPool.dblArrNCoreCols, denseCol);

         int convexity_index = m_modelCore->nBaseRows;
         //shift all the cuts over one!
         for(int r =  m_modelCore->getNumRows(); r > convexity_index; r--){
            denseCol[r] = denseCol[r-1];
         }
         denseCol[convexity_index] = 1.0; 
         //for(int r = 0; r <= m_modelCore->getNumRows(); r++){
         //  (*m_osLog) << "\ndenseCol[r:" << r << "]: " << denseCol[r];
         //}
      
         // ---
         // --- creat a sparse column from the dense column
         // ---
         sparseCol 
            = UtilPackedVectorFromDense(m_modelCore->getNumRows() + 1,
                                        denseCol, m_app->m_param.TolZero);
      }
      //(*m_osLog) << "\nPRINT sparseCol\n";
      //UtilPrintPackedVector(*sparseCol);


      DecompWaitingCol waitingCol(*li, sparseCol);
      //(*m_osLog) << "RC for this waitingCol: " 
      //           <<  m_var->getOriginalCost() - m_col->dotProduct(u);

      //what if we have created a variable that is already in LP!
      //that should not happen the way we currently do things... if it 
      //does, it could flag that we have an error of some sort
    
      //TOOD: since DecompVarList does not have its own class...
      //this is ugly, fix this later... make a helper funciton of DecompVar?

      //TODO: this is very expensive - use hash like in cuts
      if(m_varpool.isDuplicate(m_vars, waitingCol)){
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                    (*m_osLog) << "Duplicate variable, already in vars!!";
                    );
         waitingCol.deleteVar();
         waitingCol.deleteCol();
         //THINK
         m_varsThisCall--;
         m_varsThisRound--;
         // assert(0);//this can happen if multi-poly?
         continue;
      }


      if(m_varpool.isDuplicate(waitingCol)){
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                    (*m_osLog) << "Duplicate variable, already in var pool";
                    );
         waitingCol.deleteVar();
         waitingCol.deleteCol();
         //THINK
         m_varsThisCall--;
         m_varsThisRound--;
      }
      else{      
         //didn't you already do this part when you added the var to new_vars?
         //waitingCol.setReducedCost(pi, phase);
         m_varpool.push_back(waitingCol);
      }
   }
   UTIL_DELARR(denseCol);
   //printf("\nm_varpool.size(): %d", m_varpool.size());
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- addVarsToPool()       ---->";
              );  
}

// ------------------------------------------------------------------------- //
//THINK: why are these not member functions of varpool?
//why of algo? because different for different algos...
void DecompAlgo::addVarsFromPool(){

   //TODO: we have checked to make sure there are no dups added to pool
   //do we also need to check that no dup vars are added to LP? for that
   //we'd have to check across m_vars


   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- addVarsFromPool()     ---- ";
              );

   //const int maxvars_toadd = m_app->m_param.maxvars_periter;  
   //int n_newcols = std::min<int>(m_varpool.size(), maxvars_toadd);
   int n_newcols = static_cast<int>(m_varpool.size());
   //printf("\nn_newcols: %d", n_newcols);
   if(n_newcols == 0){
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << m_classTag << "  ---- addVarsFromPool()     ---->";
                 );  
      return;
   }

   //---
   //--- sort the pool by increasing reduced cost
   //---
   //THINK: need a copy constructor for object here?
   partial_sort(m_varpool.begin(), 
                m_varpool.begin() + n_newcols,
                m_varpool.end(), 
                is_less_thanD());
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 10,
              (*m_osLog) << "\nVAR POOL BEFORE:\n";
              m_varpool.print(m_osLog);
              (*m_osLog) << "\nVARS BEFORE:\n";
              printVars(m_osLog);  
              );
  
   //THINK: find_first would be good, but we need the index of the iterator
   //how do i get that?
  
   //---
   //--- never add anything with pos rc
   //---
   int index = 0;
   DecompVarPool::iterator vi;
   DecompVarPool::iterator viLast;
   for(vi = m_varpool.begin(); vi != m_varpool.end(); vi++){
      if((*vi).getReducedCost() >  -0.0000001){//TODO
         //if((*vi).getReducedCost() >  -m_app->m_param.dualTol){
         break;
      }
      index++;
   }
   n_newcols = std::min<int>(n_newcols, index);
  
   //TODO
   /*if(n_newcols > 0)
     m_cutpool.setRowsAreValid(false);*/

   //see CoinBuild

   //---
   //--- 1.) build up the block of columns to be added to the master
   //---     create a block for speed, rather than one column at a time
   //--- 2.) copy the var pointers to the DecompModel var list
   //---
   double * clb = new double[n_newcols];
   double * cub = new double[n_newcols];
   double * obj = new double[n_newcols];
   const CoinPackedVectorBase ** colBlock =
      new const CoinPackedVectorBase*[n_newcols];  
   index = 0;
   for(vi = m_varpool.begin(); vi != m_varpool.end(); vi++){
      if(index >= n_newcols)
         break;
      //THINK: bother using CoinPackedVector??
      const CoinPackedVector * col = (*vi).getColPtr();
      assert(col);
      colBlock[index] = col;
      clb[index]      = (*vi).getLowerBound();
      cub[index]      = (*vi).getUpperBound();
      obj[index]      = (*vi).getOrigCost();
      appendVars((*vi).getVarPtr());
      index++;
   }
   viLast = vi;
   m_masterSI->addCols(n_newcols, colBlock, clb, cub, obj);  
   //m_masterSI->writeMps("afterAddCol");
 
   //---
   //--- 3.) delete the col memory and clear the var pointer from varpool
   //---     the column memory is no longer needed, it has been copied into 
   //---     the master object, the variable memory is still needed, its 
   //---     pointer is now in m_vars, and no longer is needed in varpool
   //---
   //THINK is this all neccessary? just to keep memory small? or 
   //doing this for some reason of efficiency?
   for(vi = m_varpool.begin(); vi != viLast; vi++){
      (*vi).deleteCol();    
      (*vi).clearVar(); //needed? dangling pointer if not
   }
   //TODO: is this slow for vector? if so, maybe list is still the way to go
   m_varpool.erase(m_varpool.begin(), viLast);

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 10,
              (*m_osLog) << "\nVAR POOL AFTER:\n";
              m_varpool.print(m_osLog);
              (*m_osLog) << "\nVARS AFTER:\n";
              printVars(m_osLog);  
              );

   //---
   //--- free local memory
   //---
   UTIL_DELARR(colBlock);
   UTIL_DELARR(clb);
   UTIL_DELARR(cub);
   UTIL_DELARR(obj);

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- addVarsFromPool()     ---->";
              );  
}


/*-------------------------------------------------------------------------*/
//this is very close to PC code, minus reformRow, maybe we can push this
//back to base... thINK
void DecompAlgo::addCutsToPool(const double  *  x,
                               DecompCutList & newCuts,
                               int           & m_cutsThisCall){
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- addCutsToPool()       ---- ";
              );

   //for RC, make sure no cuts we are about to add are already in m_modelCore
   //also check that we have no duplicate cuts being put in here

   //TODO: do something similiar to check for pos-rc vars
   unsigned int i;
   int  r, cutIndex = 0;
   bool isViolated = false; 
   bool isDupCore;//also check relax?
   bool isDupPool;
   bool addCut;

   DecompCutPool::iterator ci;
   DecompCutList::iterator li = newCuts.begin();
   while(li != newCuts.end()){
      CoinPackedVector * row       = new CoinPackedVector();
        
      //---
      //--- create a row (in terms of original formulation, x), from a cut
      //---
      (*li)->expandCutToRow(row);

      //---
      //--- set the hash string (for quick duplicate checks)
      //---
      (*li)->setStringHash(row);
        
      bool isOptViolated = false;
      for(i = 0; i < m_optPoint.size(); i++){
         isOptViolated = (*li)->calcViolation(row, &m_optPoint[i][0]);
         if(isOptViolated)
            (*m_osLog) << "\n\nCUT VIOLATES OPT POINT";
         (*li)->print();
         assert(!isOptViolated);
      }

      //here we will use a hash table - or just STL map
      addCut    = true;
   
      isDupCore = false;
      for(r = 0; r < m_modelCore->getNumRows(); r++){
         //override isSame( )
         //in one case you check hash if expanded
         //in user case you check isSame directly
         //this will become hash lookup code
         if(m_modelCore->rowHash[r] == (*li)->getStrHash()){
            (*m_osLog) << "\n\nCUT IS DUPLICATE with Core";
            (*li)->print();
            isDupCore = true;
            break;
         }
      }

      if(isDupCore)
         addCut = false;
      else{    

         //---
         //--- is this cut already in pool 
         //---
         isDupPool = false;
         for(ci = m_cutpool.begin(); ci != m_cutpool.end(); ci++){
            if((*li)->getStrHash() == (*ci).getCutPtr()->getStrHash()){
               (*m_osLog) << "\n\nCUT IS DUPLICATE with Pool";
               (*li)->print();
               isDupPool = true;
               break;
            }           
         }

         if(isDupPool)
            addCut = false;
         else{

            isViolated = (*li)->calcViolation(row, x);//also sets it
            if(!isViolated)
               addCut = false;
         }
      }
      
      if(addCut){
         DecompWaitingRow waitingRow(*li, row);
         m_cutpool.push_back(waitingRow);
         li++;
      }
      else{
         //---
         //--- cut is not violated, do not put in cut pool, delete memory
         //---      
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                    (*m_osLog) << "\nCUT " << cutIndex 
                    << " do not put in pool";
                    );                 
         UTIL_DELPTR(*li); //need to do?
         li = newCuts.erase(li); //does this call cut destructor?
         m_cutsThisCall--;
         //then don't increment li next iter?? think
      }
      cutIndex++;
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- addCutsToPool()       ---->";
              );
   assert(m_cutsThisCall >= 0);
}



/*--------------------------------------------------------------------------*/
int DecompAlgo::addCutsFromPool(){

   //this is exactly the same for RC, except that we also have to add 
   //multipliers 

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- addCutsFromPool()     ---- ";
              );

   //TODO: do some work here to check for duplicate cuts (actually do that
   //in addCutsToPool) in RC, can add cuts that already did (no "core model"
   //)

   //TODO: do some work here to check for parallel cuts

   //TODO partial sort!
   sort(m_cutpool.begin(),        
        m_cutpool.end(), 
        is_greater_thanD());
  
#if 0
   if(m_app->m_param.logLevel > 10){
      cout << "\nCUT POOL BEFORE:" << endl;
      m_cutpool.print();
    
      //write a func for this, member of DecompModel
      cout << "\nMODEL CUTS BEFORE: " << endl;
      DecompCutList::iterator it;
      int row_index = 0;
      for(it = m_model.cuts.begin(); it != m_model.cuts.end(); it++){
         cout << row_index << " : ";
         (*it)->print();
         cout << endl;
         row_index++;
      }
      cout << endl;
   }
#endif

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nCUT POOL BEFORE:\n";
              m_cutpool.print(m_osLog);
              (*m_osLog) << "\nCUTS BEFORE:\n";
              printCuts(m_osLog);  
              );


   const int maxcuts_toadd = 100;//m_app->m_param.cut_maxcuts_periter;
  
   int n_newrows = CoinMin(static_cast<int>(m_cutpool.size()), maxcuts_toadd);
   //since we use a list - find_first won't help as it returns an 
   //iterator not an index in the list... UGH
   int index = 0;
   DecompCutPool::iterator li;
   for(li = m_cutpool.begin(); li != m_cutpool.end(); li++){
      if((*li).getViolation() < 0.00001){//PARM
         break;
      }
      index++;
   }
   //never add anything not violated
   n_newrows = std::min<int>(n_newrows, index);

   //if(n_newrows > 0)
   //  m_varpool.setColsAreValid(false);

   //TODO: look into coin build...
   double * rlb = new double[n_newrows];
   double * rub = new double[n_newrows];
   //const CoinPackedVectorBase ** rowReformBlock =
   //  new const CoinPackedVectorBase*[n_newrows];
   const CoinPackedVectorBase ** rowBlock   =
      new const CoinPackedVectorBase*[n_newrows];
  
   index = 0;
   for(li = m_cutpool.begin(); li != m_cutpool.end(); li++){
      if(index >= n_newrows)
         break;
    
      //const CoinPackedVector * rowReform = (*li).getRowReformPtr();
      const CoinPackedVector * row       = (*li).getRowPtr();
      //rowReformBlock[index] = rowReform;
      rowBlock[index]       = row;
      rlb[index] = (*li).getLowerBound();
      rub[index] = (*li).getUpperBound();
      m_cuts.push_back((*li).getCutPtr());
      m_modelCore->rowHash.push_back((*li).getCutPtr()->getStrHash());
      index++;
   }

   //actually - if we are "tightning"... we want to add to core as well!
   //THINK...
   if((m_algo == RELAX_AND_CUT) || m_isTightenAlgo) {
      //this is what we want to update in RC, but not in C
      m_modelCore->M->appendRows(n_newrows, rowBlock);
    
      //create a m_modelCore->appendRows does all of this direct from 
      //a cut pool
    
      //THINK: francois idea, when add cuts to P'?
      char   sense;
      double rhs, range;
      for(index = 0; index < n_newrows; index++){
         printf("\nadding to rowLB: %g", rlb[index]);
         printf("\nadding to rowUB: %g", rub[index]);
         m_modelCore->rowLB.push_back(rlb[index]);
         m_modelCore->rowUB.push_back(rub[index]);

         UtilBoundToSense(rlb[index], rub[index], DecompInf, 
                          sense, rhs, range);
         m_modelCore->rowRhs.push_back(rhs);
         m_modelCore->rowSense.push_back(sense);     
      }
   }
   else{
      m_masterSI->addRows(n_newrows, rowBlock, rlb, rub);
   }
   if(m_isTightenAlgo)
      m_masterSI->addRows(n_newrows, rowBlock, rlb, rub);
   //any reason to update this copy?? THINK
   //cuts on cuts... 

  
   //clean up
   index = 0;
   for(li = m_cutpool.begin(); li != m_cutpool.end(); li++){
      if(index >= n_newrows)
         break;
      //(*li).deleteRowReform();
      (*li).deleteRow();
      (*li).clearCut();//need to do this?
      index++;
   }
   m_cutpool.erase(m_cutpool.begin(), li);
  
   //UTIL_DELARR(rowReformBlock);
   UTIL_DELARR(rowBlock);
   UTIL_DELARR(rlb);
   UTIL_DELARR(rub);

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nCUT POOL AFTER:\n";
              m_cutpool.print(m_osLog);
              (*m_osLog) << "\nCUTS AFTER:\n";
              printCuts(m_osLog);  //??? add to cuts?
              (*m_osLog) << "n_newrows = " << n_newrows;
              );

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- addCutsFromPool()     ---->";
              );

   return n_newrows;
}

//TODO: compressRows, Cols


// ------------------------------------------------------------------------- //
DecompStat DecompAlgo::solutionUpdate(const DecompPhase phase,
                                      const int         maxInnerIter,
                                      const int         maxOuterIter){

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- solutionUpdate() ----\n ";
              );

   m_stats.timerOther1.restart();

   //
   // --- TODO: take maxIter steps of simplex, specific to CLP?
   // --- can you set this in CPX? we are definitely going to want this
   // --- to work with CPX sometime soon too
   //
   DecompStat stat = STAT_UNKNOWN;
  
   m_masterSI->setIntParam(OsiMaxNumIteration, maxInnerIter);
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,             
              (*m_osLog) << "Going into master solve n_cols: "
              << m_masterSI->getNumCols() << " n_rows: " 
              << m_masterSI->getNumRows() << endl;
              );

   //TODO?
   //THINK: if last iteration was price, do interior?
   //if last iteration was cut, do simplex? think about how
   //basis gets updated, etc.. hmmm... could be interesting

  
   //if we allow for interior, need crossover too?
#ifdef __DECOMP_LP_CPX__
   CPXENVptr env
      = dynamic_cast<OsiCpxSolverInterface*>(m_masterSI)->getEnvironmentPtr();
   CPXLPptr  lp
      = dynamic_cast<OsiCpxSolverInterface*>(m_masterSI)->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
   CPXsetintparam( env, CPX_PARAM_PREIND, CPX_ON );
   CPXsetintparam( env, CPX_PARAM_SCRIND, CPX_ON );
   CPXsetintparam( env, CPX_PARAM_SIMDISPLAY, 2 );
#endif
  
   switch(phase){
   case PHASE_INIT:
      //
      // --- solve with primal simplex? or dual??
      //    
      if(m_algo == DECOMP)
         m_masterSI->setHintParam(OsiDoPresolveInInitial, false, OsiHintDo);
#ifdef DO_INTERIOR
      CPXhybbaropt(env, lp, 0);
#else
      if(m_nodeIndex == 0)
         m_masterSI->initialSolve();
      else
         m_masterSI->resolve();
#endif
      break;
   case PHASE_PRICE:
      //m_si->setHintParam(OsiDoDualInResolve, false);
      //otherwise hit dual cutoff? default solves to optimal, i think
      //this is only useful for in branch and bound
    
      m_masterSI->setDblParam(OsiDualObjectiveLimit, DecompInf);

      if(m_algo == DECOMP)//THINK!
         m_masterSI->setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
    
      //m_masterSI->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
      //force primal??? probably yes...  THINK about this
#ifdef DO_INTERIOR
      CPXhybbaropt(env, lp, 0);
#else
      //can't use primal if stats_infeasible, need dual ray...
      //STOP - try using primal - reduce number of simplex iters
      //if(!m_masterSI->isProvenPrimalInfeasible())
      //   m_masterSI->setHintParam(OsiDoDualInResolve, false, OsiHintDo);
      m_masterSI->resolve();//THINK: want to use dual or primal
#endif
      break;
   case PHASE_CUT:
#ifdef DO_INTERIOR
      CPXhybbaropt(env, lp, 0);
#else
      m_masterSI->resolve();//THINK: want to use dual or primal
#endif
      break;
   default:
      assert(0);
   }
  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << "\nisAbandoned()                 : " 
              << m_masterSI->isAbandoned();
              (*m_osLog) << "\nisProvenOptimal()             : " 
              << m_masterSI->isProvenOptimal();
              (*m_osLog) << "\nisProvenPrimalInfeasible()    : " 
              << m_masterSI->isProvenPrimalInfeasible();
              (*m_osLog) << "\nisProvenDualInfeasible()      : " 
              << m_masterSI->isProvenDualInfeasible();
              (*m_osLog) << "\nisPrimalObjectiveLimitReached : " 
              << m_masterSI->isDualObjectiveLimitReached();
              (*m_osLog) << "\nisDualObjectiveLimitReached   : " 
              << m_masterSI->isDualObjectiveLimitReached();
              (*m_osLog) << "\nisIterationLimitReached       : " 
              << m_masterSI->isIterationLimitReached();
              );
             
             
   //THINK - about all possible return codes, CLP second stat needed?
   //what if stopped on simplex iterations?

   //none of these things are updated if you use interior - you have to
   //do it manually... THINK - design, this section shouldn't depend on
   //OSI or any solver - it should just ask for primal and dual solution
   //info - and then depending on the method as param or derived, it
   //goes and gets what it needs... THINK

   //STOP interior
   if(m_masterSI->isProvenOptimal()){
      stat = STAT_FEASIBLE;  
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,             
                 const double * masterSol = m_masterSI->getColSolution();
                 for(int i = 0; i < m_masterSI->getNumCols(); i++){
                    if(fabs(masterSol[i]) > m_app->m_param.TolZero) 
                       (*m_osLog) << "\nmasterSol[" << i << "]: " << masterSol[i];
                 }
                 /*
                   const double * act   = m_masterSI->getRowActivity();
                   const double * rhs   = m_masterSI->getRightHandSide();
                   const char   * sense = m_masterSI->getRowSense();
                   for(int i = 0; i < m_masterSI->getNumRows(); i++){
                   (*m_osLog) << "\nrow: " << i << " " << act[i] << " "
                   << sense[i] << " " << rhs[i];                   
                   }
                 */
                 );      

      //sanity check: cx = ub + alpha
      double primalObj = m_masterSI->getObjValue();
      double dualObj   = 0.0;
      const double * rhs  = m_masterSI->getRightHandSide();
      const double * dual = m_masterSI->getRowPrice();
      for(int i = 0; i < m_masterSI->getNumRows(); i++){
         dualObj += dual[i] * rhs[i];
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,             
                    printf("\ni: %d, dualObj: %g, dual: %g, rhs: %g",
                           i, dualObj, dual[i], rhs[i]);
                    );
      }
      //Dmitry notes...
      //---
      /*
        for (i = 0; i < rows; i++) {
        dualObj += y[i] * b[i];
        }

        for (i = 0; i < columns; i++) {
        ch = boundType[i];

        if (ch & BOUND_LOWER) {
        dualObj += d1[i] * l[i];
        }

        if (ch & BOUND_UPPER) {
        dualObj -= d2[i] * u[i];
        }
        }
      */
     
      //--- BUT, for our master, we are doing
      //--- 0 <= lambda <= 1... which means we need d2, but an LP
      //--- solver never provides this... but we do NOT need that
      //--- ub on lambda!! since we have the convexity constraint!
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,             
                 (*m_osLog) << "\nPrimal obj: " << primalObj
                 << " Dual obj: " << dualObj;
                 );
      //should this be expected??
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,             
                 if(!UtilIsZero(primalObj-dualObj, 1.0e-1)){       
                    (*m_osLog) << "OFF BY" << fabs(primalObj - dualObj);
                 }
                 );
      //assert(UtilIsZero(primalObj-dualObj, 1.0e-1));
    
   }
   else if(m_masterSI->isProvenPrimalInfeasible()){
      stat = STAT_INFEASIBLE;

      //if this is the case, we have to get a dual ray - but if presolve
      //determines it infeasible, we need to go in and force simplex to do so
      //any way to tell if it was presolve or simplex that found in inf?
      //from status code? CPXERR_PRESLV_INForUNBD?
      m_masterSI->setHintParam(OsiDoPresolveInResolve, false, OsiHintDo);
      m_masterSI->resolve();
      m_masterSI->setHintParam(OsiDoPresolveInResolve, true, OsiHintDo);

    
   }
   else{
#ifdef __DECOMP_LP_CLP__
      printf("\nCLP stat: %d, 2ndStat: %d", 
             m_masterCLP->problemStatus(), 
             m_masterCLP->secondaryStatus()
             );
#endif
      assert(0);
   }

   
   //TODO: add to stats
   //how many simplex iterations? are we warm-starting ok?
   //printf("\nSolUpdate ITER: %d", m_masterSI->getIterationCount());

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- solutionUpdate() ----> ";
              );

   m_stats.thisSolUpdate.push_back(m_stats.timerOther1.timeElapsed()); 
  
   return stat;
}

// ------------------------------------------------------------------------ //
DecompPhase DecompAlgo::phaseUpdate(const DecompPhase phase,
                                    const DecompStat  stat){
  
   bool        isCutPossible, isPricePossible, mustSwitch, considerSwitch;
   DecompPhase nextPhase = PHASE_UNKNOWN;
 
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- phaseUpdate() ---- ";
              );  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,  
              (*m_osLog) << "\nm_cutsThisRound  : " << m_cutsThisRound;
              (*m_osLog) << "\nm_varsThisRound  : " << m_varsThisRound;
              (*m_osLog) << "\nm_cutsThisCall   : " << m_cutsThisCall;
              (*m_osLog) << "\nm_varsThisCall   : " << m_varsThisCall;
              (*m_osLog) << "\nm_cutCallsTotal  : " << m_cutCallsTotal;
              (*m_osLog) << "\nm_priceCallsTotal: " << m_priceCallsTotal;
              (*m_osLog) << "\nm_cutCallsRound  : " << m_cutCallsRound;
              (*m_osLog) << "\nm_priceCallsRound: " << m_priceCallsRound;
              (*m_osLog) << "\nPHASEIN          : " << DecompPhaseStr[phase];
              (*m_osLog) << "\nSTATIN           : " << DecompStatStr[stat];
	      (*m_osLog) << "\nm_tlb            : " << m_tlb;
	      (*m_osLog) << "\nm_tub            : " << m_tub;
              );

   //---
   //--- if the true lower bound meets the global ub, we are done
   //---
   //m_tlb might not be calculated at this point
   //if(m_algo == PRICE_AND_CUT && phase != PHASE_INIT){
      if(m_tlb >= (m_tub - 1.0e-5)){
	 nextPhase = PHASE_DONE;
	 goto PHASE_UPDATE_FINISH;
      }
      //}

   //new
   //app user feasible does not need to be virtual
   //the default checks against core and relaxed constraints - see SmallIP
   bool ipFeasible;
   bool appFeasible;
   if(ipFeasible = isIPFeasible(m_xhat)){
      //printf("\nIP FEASIBLE, we are done? there might still be cuts"); 
      if(appFeasible = m_app->APPisUserFeasible(m_xhat,
						m_modelCore->getNumCols(),
						m_app->m_param.TolZero)){
	 //printf("\nAPP FEASIBLE"); 
      }        
      if(ipFeasible && !appFeasible){
	 nextPhase = PHASE_CUT;  
	 goto PHASE_UPDATE_FINISH;
      }
   }
   
   //---
   //--- we have exceeded the cut iter limit and the price iter limit
   //--- we are done
   //---
   if((m_priceCallsTotal >= m_app->m_param.LimitTotalPriceIters) &&
      (m_cutCallsTotal   >= m_app->m_param.LimitTotalCutIters)){
      nextPhase = PHASE_DONE;  
      goto PHASE_UPDATE_FINISH;
   }

   //what if user puts in total=1, round=0?
   isCutPossible   = (m_cutCallsTotal   < m_app->m_param.LimitTotalCutIters);
   isPricePossible = (m_priceCallsTotal < m_app->m_param.LimitTotalPriceIters);
   //printf("\nisCutPossible: %d", isCutPossible);
   //printf("\nisPricePossible: %d", isPricePossible);
  
   if(stat == STAT_INFEASIBLE){
      //what if we are not pricing and it goes infeasible - that should STOP
      //test this with an infeasible milp example
    
      //deal with the case of really infeasible
      if(phase == PHASE_CUT){
         nextPhase = PHASE_PRICE;
      } else if((m_priceCallsRound > 0) && (m_varsThisCall == 0)){
         nextPhase = PHASE_DONE;
      }
      else{
         nextPhase = PHASE_PRICE;
      }
      goto PHASE_UPDATE_FINISH;
   }
  
   //TODO: what if goes infeasible in middle?

   switch(phase){
   case PHASE_INIT:
      {
         if(isPricePossible)
            nextPhase = PHASE_PRICE;
         else if(isCutPossible)
            nextPhase = PHASE_CUT;
         else
            nextPhase = PHASE_DONE;
      }
      break;
   case PHASE_PRICE:
      {
         mustSwitch     = false;
         considerSwitch = false;    
      
         if(!isPricePossible || (m_varsThisCall == 0) || (m_varsThisRound == 0))
            mustSwitch = true;
         if(m_priceCallsRound >= m_app->m_param.LimitRoundPriceIters)
            considerSwitch = true;
      
         if(mustSwitch){
            //---
            //--- we must switch from pricing
            //---
            if(!isCutPossible){
               //---
               //--- if we exceed the cut iter limit, we are done
               //---
               nextPhase = PHASE_DONE;
            }
            else{
               if((m_cutCallsTotal > 0)  && 
                  (m_cutsThisRound == 0) && 
                  (m_varsThisRound == 0)){
                  //---
                  //--- nothing new happened, so we are done
                  //---
                  nextPhase = PHASE_DONE;
               }
               else{
                  //---
                  //--- something new happened, so try cuts again
                  //---
                  nextPhase = PHASE_CUT;
                  m_cutsThisRound = 0;
                  m_cutCallsRound = 0;
               }
            }
         }//END: if(mustSwitch)
         else if(considerSwitch){
            //---
            //--- we consider switching from pricing
            //---
            if(!isCutPossible){
               if(!isPricePossible){
                  //---
                  //--- if we exceed both iter limits, we are done
                  //---
                  nextPhase = PHASE_DONE;
               }
               else{
                  //---
                  //--- if we exceed the cut iter limit, but not the price limit
                  //--- since we are not in mustSwitch, m_varsThisRound > 0, so 
                  //--- we can go back to pricing, even though it violates the 
                  //--- round counter, because we have no other choice
                  //---
                  nextPhase = PHASE_PRICE;
               }
            }
            else{
               if((m_cutsThisRound == 0) && (m_varsThisRound == 0)){
                  //---
                  //--- nothing new happened, so we are done
                  //---
                  nextPhase = PHASE_DONE;
               }
               else{
                  //---
                  //--- something new happened, so try cuts again
                  //---
                  nextPhase = PHASE_CUT;
                  m_cutsThisRound = 0;
                  m_cutCallsRound = 0;
               }
            }
         } //END: else if(considerSwitch)
         else
            nextPhase = PHASE_PRICE;
      }

      break;
   case PHASE_CUT:
      {
         mustSwitch     = false;
         considerSwitch = false;    
      
         if(!isCutPossible || (m_cutsThisCall == 0) || (m_cutsThisRound == 0))
            mustSwitch = true;
         if(m_cutCallsRound >= m_app->m_param.LimitRoundCutIters)
            considerSwitch = true;
      
         if(mustSwitch){
            //---
            //--- we must switch from cutting
            //---
            if(!isPricePossible){
               //---
               //--- if we exceed the price iter limit, we are done
               //---
               nextPhase = PHASE_DONE;
            }
            else{
               if((m_cutsThisRound == 0) && (m_varsThisRound == 0)){
                  //---
                  //--- nothing new happened, so we are done
                  //---
                  nextPhase = PHASE_DONE;
               }
               else{
                  //---
                  //--- something new happened, so try price again
                  //---
                  nextPhase = PHASE_PRICE;
                  m_varsThisRound   = 0;
                  m_priceCallsRound = 0;
               }
            }
         }//END: if(mustSwitch)
         else if(considerSwitch){
            //---
            //--- we consider switching from cutting
            //---
            if(!isPricePossible){
               if(!isCutPossible){
                  //---
                  //--- if we exceed both iter limits, we are done
                  //---
                  nextPhase = PHASE_DONE;
               }
               else{
                  //---
                  //--- if we exceed the price iter limit, but not the cut limit
                  //--- since we are not in mustSwitch, m_cutsThisRound > 0, so 
                  //--- we can go back to cutting, even though it violates the 
                  //--- round counter, because we have no other choice
                  //---
                  nextPhase = PHASE_CUT;
               }
            }
            else{
               if((m_cutsThisRound == 0) && (m_varsThisRound == 0)){
                  //---
                  //--- nothing new happened, so we are done
                  //---
                  nextPhase = PHASE_DONE;
               }
               else{
                  //---
                  //--- something new happened, so try price again
                  //---
                  nextPhase = PHASE_PRICE;
                  m_varsThisRound   = 0;
                  m_priceCallsRound = 0;
               }
            }
         } //END: else if(considerSwitch)
         else
            nextPhase = PHASE_CUT;
      }
      break;
   case PHASE_DONE:
      break;
   case PHASE_UNKNOWN:
   default:
      assert(0);//TODO: real error codes?
   }

 PHASE_UPDATE_FINISH:
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,  
              (*m_osLog) << "\nPHASEOUT         : "
              << DecompPhaseStr[nextPhase] << "\t";
              );   
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- phaseUpdate() ----> ";
              );
   return nextPhase;
}

#if 0
// ------------------------------------------------------------------------ //
DecompPhase DecompAlgo::phaseUpdate(const DecompPhase   phase,
                                    const DecompStat    stat,
                                    const int           n_cutCallsTotal,
                                    const int           n_priceCallsTotal,
                                    int               & m_cutsThisCall,
                                    int               & m_varsThisCall,
                                    int               & n_cutCalls,
                                    int               & n_priceCalls){
  
   //why not part of DecompAlgo object? n_... 

   //THINK: rather than call it price,cut - call it inner outer?

   //TODO: need a frequency parameter for how you do things...
   //a) which to start with, b) how many consecutive price, cut
   //or allow it to price out, then cut till none... many strategies to offer
   //need tailing off detection

   //options:
   // just price = (a), x = inf, y = 0
   // just cut   = (a), x = 0,   y = inf
   // if that is the case, then phaseUpdate in PC can be base... 
   // just set the correct (a) and (b)
   // (a) price out up to x iters, then cut up to y iters, repeat
   //   for now just look at (a)
   // (b) cut out up to x iters, then price up to y iters

   //---
   //--- Consider several different scenarios.
   //---
   //--- (1)
   //---  cutIterRound   = 1
   //---  priceIterRound = 1
   //---  cutIterTotal   = 5
   //---  priceIterTotal = 5
   //---
   //---   init                   cutRound cutTotal priceRound priceTotal
   //---   price -> infeasible       0        0          1         1
   //---   price -> infeasible       0        0          2         2
   //---   price -> found var        0        0          3*        3
   //---   cut   -> found cut        1*       1          0         3
   //---   price -> found var        0        1          1*        4
   //---   cut   -> found cut        1*       2          0         4*
   //---   price -> found var        0        2          1*        5*
   //---   cut   -> found cut        1        3          0         5*
   //---   cut   -> found cut        2        4          0         5
   //---   cut   -> no cut found*    3        5          0         5
   //---   done
   //---
   //--- (2)
   //---  cutIterRound   = 1
   //---  priceIterRound = 1
   //---  cutIterTotal   = 5
   //---  priceIterTotal = 3
   //---
   //---   init                   cutRound cutTotal priceRound priceTotal
   //---   price -> infeasible       0        0          1         1
   //---   price -> infeasible       0        0          2         2
   //---   price -> infeasible       0        0          3         3*
   //---   done
   //---
   //--- (3)
   //---  cutIterRound   = 2
   //---  priceIterRound = 2
   //---  cutIterTotal   = 6
   //---  priceIterTotal = 8
   //---
   //---   init                   cutRound cutTotal priceRound priceTotal
   //---   price -> found var        0        0          1         1
   //---   price -> found var        0        0          2*        2
   //---   cut   -> found cut        1        1          0         2
   //---   cut   -> found cut        2*       2          0         2
   //---   price -> found var        0        2          1         3
   //---   price -> found var        0        2          2*        4
   //---   cut   -> no cut found*    1        3          0         4
   //---   price -> found var        0        3          1         5
   //---   price -> no var found*    0        3          2         6
   //---   cut   -> cut found        1        4          0         6
   //---   cut   -> no cut found*    2        5          0         6
   //---   price -> no var found*    2        5          1         7
   //---   done

   //cut   -> no cut found
   //price -> found var
   //price -> no var found
   //should cut, not be done!


   DecompPhase nextPhase = PHASE_UNKNOWN;
#if 0  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- phaseUpdate() ---- ";
              );


   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,  
              (*m_osLog) << "\nm_cutsThisCall: "   << m_cutsThisCall;
              (*m_osLog) << " m_varsThisCall: "    << m_varsThisCall;
              (*m_osLog) << " n_cutCallsTotal: "   << n_cutCallsTotal;
              (*m_osLog) << " n_priceCallsTotal: " << n_priceCallsTotal;
              (*m_osLog) << " n_cutCalls: "   << n_cutCalls;
              (*m_osLog) << " n_priceCalls: " << n_priceCalls;
              (*m_osLog) << "\nPHASEIN : "    << DecompPhaseStr[phase] << "\t";
              (*m_osLog) << "STAT IN: "       << DecompStatStr[stat];
              );

   //stopping criteria - different override for each algo

   //---
   //--- this takes care of the case where we have exceeded the total
   //--- number of iterations for both price and cut, in this case, we are done
   //---
   if((n_priceCallsTotal >= m_app->m_param.LimitTotalPriceIters)
      && (n_cutCallsTotal >= m_app->m_param.LimitTotalCutIters)){
      nextPhase = PHASE_DONE;
      goto PHASE_UPDATE_FINISH;
   }
    
   switch(phase){
   case PHASE_INIT:
      //THINK about this
      if(m_algo == CUT)
         nextPhase = PHASE_CUT;      
      else if(n_priceCalls < m_app->m_param.LimitRoundPriceIters)
         nextPhase = PHASE_PRICE;
      else if(n_cutCalls < m_app->m_param.LimitRoundCutIters)
         nextPhase = PHASE_CUT;
      else
         nextPhase = PHASE_DONE;
      break;
  
   case PHASE_PRICE:
      //---
      //--- we are pricing
      //---
      if(stat == STAT_FEASIBLE){
         //---
         //--- we have a feasible solution 
         //---
         if(n_priceCallsTotal < m_app->m_param.LimitTotalPriceIters){
            if(n_priceCalls < m_app->m_param.LimitRoundPriceIters){
               if(m_varsThisCall > 0){
                  //---
                  //--- we have not exceeded either price iter limit and
                  //--- we have not priced out, continue pricing
                  //---
                  nextPhase = PHASE_PRICE;
               }           
               else{
                  //---
                  //--- we have not exceeded either price iter limit and
                  //--- but, we priced out, consider cutting or done
                  //---
                  if(n_cutCallsTotal < m_app->m_param.LimitTotalCutIters){
                     if(m_cutsThisCall > 0){
                        //---
                        //--- we have not exceeded the cut iter limit and
                        //--- we found a cut in the last cut pass, switch to cut
                        //---
                        nextPhase  = PHASE_CUT;
                        m_cutsThisCall  = 0;
                        n_cutCalls = 0;
                     }
                     else{
                        //---
                        //--- we have priced out and found no cuts last iter
                     }
                  }
               }
            }
         }
      




         if(n_priceCallsTotal < m_app->m_param.LimitTotalPriceIters){
            if(n_priceCalls < m_app->m_param.LimitRoundPriceIters){
               //---
               //--- 
               //---
               nextPhase = PHASE_PRICE;
               //but have to check if we generated any vars... if not,
               //same questions as below... make separate func?
            }
            else{
               //---
               //--- we have exceeded the price iter limit this round
               //---
               if(n_cutCallsTotal < m_app->m_param.LimitTotalCutIters){
                  if(n_cutCalls < m_app->m_param.LimitRoundCutIters){
                     //---
                     //--- we have not exceeded either cut iter limit
                     //---
                     if((n_cutCallsTotal == 0) || (m_cutsThisCall > 0)){
                        //---
                        //--- the last cut round (if exists) was successful 
                        //---   switch to cut, reset counters
                        //---
                        nextPhase  = PHASE_CUT;
                        m_cutsThisCall  = 0; 
                        n_cutCalls = 0;
                     }
                     else{
                        //---
                        //--- the last cut round was unsuccessful
                        //---   switch to price, reset counters
                        //---
                        nextPhase    = PHASE_PRICE;
                        m_varsThisCall    = 0;
                        n_priceCalls = 0;
                     }
                  }
                  else{
                     //---
                     //--- we have exceeded the cut iter limit this round
                     //---
                     nextPhase = PHASE_DONE;
                  }
               }
               else{
                  //---
                  //--- we have exceeded the cut iter limit 
                  //---
                  nextPhase = PHASE_DONE;
               }
            }
         }
         else{
            //---
            //--- we have exceed the price iter limit
            //---
            if(n_cutCallsTotal < m_app->m_param.LimitTotalCutIters){
               if(n_cutCalls < m_app->m_param.LimitRoundCutIters){
                  //---
                  //--- we have not exceeded either cut iter limit
                  //---
                  if(m_cutsThisCall > 0){
                     //---
                     //--- the last cut round was successful 
                     //---   switch to cut, reset counters
                     //---
                     nextPhase  = PHASE_CUT;
                     m_cutsThisCall  = 0; 
                     n_cutCalls = 0;
                  }
                  else{
                     //---
                     //--- the last cut round was unsuccessful
                     //---
                     if(m_varsThisCall > 0){
                        nextPhase  = PHASE_CUT;
                        m_cutsThisCall  = 0; 
                        n_cutCalls = 0;
                     }
                     else
                        nextPhase  = PHASE_DONE;
                  }
               }
               else{
                  //---
                  //--- we have exceeded the cut iter limit this round
                  //---
                  if(m_varsThisCall > 0){
                     nextPhase  = PHASE_CUT;
                     m_cutsThisCall  = 0; 
                     n_cutCalls = 0;
                  }
                  else
                     nextPhase = PHASE_DONE;
               }
               else{
                  //---
                  //--- we have exceeded the cut iter limit 
                  //---
                  nextPhase = PHASE_DONE;
               }
            }
         }
      }
      else{
         //---
         //--- master is infeasible, we have to price to break infeasiblity
         //--- unless we have hit the price iter limit, if so, we are done
         //---
         if(n_priceCallsTotal > m_app->m_param.LimitTotalPriceIters)
            nextPhase = PHASE_DONE;
         else
            nextPhase = PHASE_PRICE;
      }
      //TODO: what about problems that are really infeasible - THINK
      break;
   case PHASE_CUT:
      //---
      //--- we are cutting
      //---
      if(stat == STAT_FEASIBLE){
         if((n_cutCalls >= m_app->m_param.LimitRoundCutIters) || (m_cutsThisCall == 0)){
            //---
            //---  too many cut steps in a row OR we found no violations
            //---
            if((m_algo == CUT) || ((m_varsThisCall == 0) && (n_priceCalls > 0))){
               //---
               //--- we didn't find any vars last iteration either, DONE
               //---
               nextPhase = PHASE_DONE;
            }
            else if(n_priceCalls < m_app->m_param.LimitRoundPriceIters){            
               //--- 
               //--- switch to price phase, reset the counters
               //---
               nextPhase    = PHASE_PRICE;
               m_varsThisCall    = 0; 
               n_priceCalls = 0;
            }
            else{
               //---
               //--- we found no violations, and called too many price 
               //--- steps in a row, DONE
               //---
               nextPhase = PHASE_DONE;
            }
         }
         else{
            //---
            //--- we found a cut last time, try again
            //---
            if(m_algo == RELAX_AND_CUT){
               //after cut found, we must price next!
               //THINK...  cutIters=1
               //--- 
               //--- switch to price phase, reset the counters
               //---
               nextPhase    = PHASE_PRICE;
               m_varsThisCall    = 0; 
               n_priceCalls = 0;
            }
            else{
               nextPhase = PHASE_CUT;
            }
         }
      }
      else{
         //---
         //--- master is infeasible, we have to price to break infeasiblity
         //---
         nextPhase = PHASE_PRICE;
      }
      //TODO: what about problems that are really infeasible - THINK
      break;
   case PHASE_DONE:
      break;
   case PHASE_UNKNOWN:
   default:
      assert(0);//TODO: real error codes?
   }

 PHASE_UPDATE_FINISH:
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,  
              (*m_osLog) << "\nPHASEOUT: " 
              << DecompPhaseStr[nextPhase] << "\t";
              );   
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- phaseUpdate() ----> ";
              );
#endif
   return nextPhase;
}
#endif

#ifdef __DECOMP_LP_CLP__
#if 1
// ------------------------------------------------------------------------- //
vector<double*> DecompAlgo::getDualRays(int maxNumRays){

   vector<double*> rays = m_masterSI->getDualRays(maxNumRays);

   const double * ray = rays[0];
   CoinAssert(ray);

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,  
              bool isProof = isDualRayInfProof(ray,
                                               m_masterSI->getMatrixByRow(),
                                               m_masterSI->getColLower(),
                                               m_masterSI->getColUpper(),
                                               m_masterSI->getRightHandSide(),
                                               m_osLog);
              printBasisInfo(m_masterSI, m_osLog);
              );
   CoinAssert(isDualRayInfProof(ray,
                                m_masterSI->getMatrixByRow(),
                                m_masterSI->getColLower(),
                                m_masterSI->getColUpper(),
                                m_masterSI->getRightHandSide()));
   
   
#if 0
   const int m = m_masterSI->getNumRows();
   double * pneg = new double[m];
   transform(ray, ray + m, pneg, negate<double>());

   vector<double*> rays2;
   rays2.push_back(pneg);
   return rays2;
#else
   return rays;
#endif


#if 0
     
   int    * binds     = new int[m];    //basic indices (negative if artifical)
   double * bels      = new double[m]; //basic variable values  
   double * binvAR    = new double[n]; //row of the tableau B-1 A
   double * binvR     = new double[m]; //row of the basis inverse

   CPXgetbhead(siCpx->getEnvironmentPtr(), siCpx->getLpPtr(), binds, bels);


#ifdef DBG_PRINT
   for(int r = 0; r < m; r++){
      if(fabs(bels[r]) > param.tol)
         (*log) << "Basic Var : " << binds[r] << " Value -> " << bels[r] << endl;
   }
#endif


   //for DECOMP, colLB is all 0, so don't really need to do this... 
   const double * colLB = si->getColLower();
   int multiplier, i;
 
   for (int r = 0; r < m; r++){

      if(rays.size() >= maxNumRays) break;
    
      bool isArtificial = binds[r] < 0 ? true : false;
      if((!isArtificial && bels[r] < colLB[binds[r]] - param.tol)
         || (isArtificial && bels[r] > param.tol)){
         //if basic variable is below its lower bound 
         //OR basic variables is articifial and positive
      
         multiplier = bels[r] < -1 * param.tol ? 1 : -1;
         //can this just check isArtificial?

#ifdef DBG_PRINT
         (*log) << "row " << r << " : bels[r] : " << bels[r] 
                << " mult -> " << multiplier << endl;
#endif

         CPXbinvarow(siCpx->getEnvironmentPtr(), siCpx->getLpPtr(), r, binvAR);

#ifdef DBG_PRINT
         (*log) << "binvA[" << r << "] : ";
         for(i = 0; i < n; i++){
            if(multiplier * binvAR[i] < -1 * param.tol)
               (*log) << "-";
            else (*log) << "+";
         }
         (*log) << endl;
#endif

         for (i = 0; i < n; i++){
            //if break, then we cannot use this row to prove unboundeness:
            if (multiplier * binvAR[i] < -1 * param.tol)
               break;    
         }
         if (i < n) continue;

         //If we get here, we have found a dual ray which proves 
         //primal infeasible / dual unbounded.
      
         //Case 1: basic is below its lower bound but all the elements 
         //in the tableau are positive
         //Case 2: the artifcial variables is positive but all the elements
         //in the tableau are negative

#ifdef DBG_PRINT
         (*log) << "row " << r << " provides a proof" << endl;
#endif
    
         CPXbinvrow(siCpx->getEnvironmentPtr(), siCpx->getLpPtr(), r, binvR);
      
         double * ray = new double[m];
         for(i = 0; i < m; i++)
            ray[i] = multiplier * binvR[i];
         rays.push_back(ray);
      }
   }

   delete [] binds;
   delete [] bels;
   delete [] binvAR;
   delete [] binvR;

   return rays;
#endif
}

#else

// ------------------------------------------------------------------------- //
vector<double*> DecompAlgo::getDualRays(int maxNumRays){
   //only good for one currently

   vector<double*> rays = m_masterSI->getDualRays(maxNumRays);
   const double * ray = rays[0];
   bool isProof = isDualRayInfProof(ray,
                                    m_masterSI->getMatrixByRow(),
                                    m_masterSI->getColLower(),
                                    m_masterSI->getColUpper(),
                                    m_masterSI->getRightHandSide(),
                                    m_osLog);
   //assert(isProof);
   return rays;
}
#endif
#endif

#ifdef __DECOMP_LP_CPX__
// ------------------------------------------------------------------------- //
//TEST THIS
vector<double*> DecompAlgo::getDualRays(int maxNumRays){
   OsiCpxSolverInterface * siCpx
      = dynamic_cast<OsiCpxSolverInterface*>(m_masterSI);
   const int m = m_masterSI->getNumRows();
   const int n = m_masterSI->getNumCols();
   double proof_p;
   bool   isProof;
  
   vector<double*> rays;
   double * ray = new double[m];
   int err = CPXdualfarkas(siCpx->getEnvironmentPtr(),
                           siCpx->getLpPtr(), ray, &proof_p);//proof_p
   if(err){
      cerr << "CPXdualfarkas returns err " << err << endl;
      abort();
   }
   cout << "\nAfter dual farkas proof_p = " << proof_p;


   
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,  
              bool isProof = isDualRayInfProof(ray,
                                               m_masterSI->getMatrixByRow(),
                                               m_masterSI->getColLower(),
                                               m_masterSI->getColUpper(),
                                               m_masterSI->getRightHandSide(),
                                               m_osLog);
              printBasisInfo(m_masterSI, m_osLog);
              );
  


  
   //cpx is flip of clp?
   double * pneg = new double[m];
   transform(ray, ray + m, pneg, negate<double>());
  
   rays.push_back(pneg);
   return rays;
}
#endif


#if 1
// ------------------------------------------------------------------------- //
//THINK: PC specific? want for RC too
int DecompAlgo::generateInitVars(DecompVarList & initVars){

   unsigned int attempts;
   int         c, m_varsThisCall;
   double      aveC;
  
   const int   n_corecols = m_modelCore->getNumCols();
   double *    objCoeff   = m_app->m_model.objCoeff;

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- generateInitVars() ---- ";
              );
  
   //---
   //--- APP: create an initial set of points F'[0] subseteq F'
   //--- The base implementation of this function does nothing.
   //--- This is the user's chance to implement something application 
   //--- specific.
   //---
   m_varsThisCall = m_app->generateInitVars(initVars, m_whichModel);
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
              (*m_osLog) << "\nm_varsThisCall from app = " 
              << m_varsThisCall << "\n";
              );

   //THINK: is this somehow specific to PC? doesn't this belong in 
   //base APP implementatin of generateInitVars

   //---
   //--- APP: implementation of generatInitVars
   //--- Note, we check for duplicates here, we stop once we have
   //--- hit LimitInitVars vars, or until we've tried 2 * LimitInitVars 
   //--- times.
   //---


   //TODO: if we called tighten before - or even if not, solve the LP
   //over original formulation to generate initial points uisng last
   //optimal dual vector... also do with cost vector c, since we are
   //solving heuristics of original problem  - this could get lucky and
   //get the optimal column right off the bat
  
   if(initVars.size() < m_app->m_param.LimitInitVars){
    
      // ---
      // --- create an initial set of points F'[0] subseteq F'
      // --- randomly by solving zSP(c + eps), eps = U[0,ave(c)]
      // ---
      //create a first one with original c, think about case of using TSP-heur
      double * costeps = new double[n_corecols];
      assert(objCoeff);
      aveC = UtilAve(objCoeff, n_corecols);

      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
                 (*m_osLog) << "\naveC = " << aveC;
                 );
     
      attempts = 0;
      while((initVars.size() < m_app->m_param.LimitInitVars)
            && (attempts < (2 * m_app->m_param.LimitInitVars))){
      
         //---
         //--- perturb the cost vector
         //---
         srand(attempts);
         for(c = 0; c < n_corecols; c++){
            double r = 0.0;
            if(attempts != 0)
               r = UtilURand(-aveC, aveC);
            costeps[c] = objCoeff[c] + r;
            //printf("\ncosteps[%d]: %g -> obj: %g + rand: %g",
            //       c, costeps[c], objCoeff[c], r);
         }


         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 50,
                    string fileName = "subProbInit.a" + UtilIntToStr(attempts);
                    printCurrentProblem(m_subprobSI[m_whichModel], fileName);
                    );
                 
         //TODO: multi-poly
      
         //---
         //--- APP: solve zSP(c + eps) 
         //---
         //TODO: stat return?
         solveRelaxed(m_whichModel,
                      costeps,               //reduced cost (fake here)
                      objCoeff,              //original cost vector 
                      0.0,                   //alpha        (fake here)
                      n_corecols,            //num core columns
                      false,                 //check for rc < 0 ?
                      true,                  //check for duplicate vars
                      m_subprobSI[m_whichModel],  //subprob solver interface
                      initVars);             //var list to populate
      
         // THINK: check for duplicate variables - done in solveRelaxed
         // don't assume the user does the duplicate check - should be done
         // by col pool also
         attempts++;
      }
    
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 4,
                 (*m_osLog) << "\nm_varsThisCall = " << initVars.size() << "\n";
                 );


      //--- 
      //--- ?? solve a few iterations of subgradient to get init vars?
      //---
    
      //--- put them in the var pool    

      UTIL_DELARR(costeps);
   }

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- generateInitVars() ----> ";
              );  
   return static_cast<int>(initVars.size());
}
#endif




//TODO: solveRelaxed with multiPoly -
//perfect for parallel - future extensions? if works well.. 

// --------------------------------------------------------------------- //
//TODO: also need to accept fixed variables (for consistency
//and for branching), its possible this function can return more
//than one var
DecompStat DecompAlgo::solveRelaxed(const int             whichModel,
                                    const double        * redCostX,
                                    const double        * origCost,
                                    const double          alpha,
                                    const int             n_origCols,
                                    const bool            checkRC,
                                    const bool            checkDup,
                                    OsiSolverInterface  * subprobSI,
                                    list<DecompVar*>    & vars){

   // ---
   // --- For pricing,
   // --- redCostX: is the red-cost for each original column  (c - uhat A")_e
   // --- origCost: is the original cost for each original column c_e
   // --- alpha:    is the dual for the convexity constraint 
   // --- 
   // --- The reduced cost of a new variable (column) is the sum of the 
   // --- reduced cost on each of the original columns in the new variable
   // --- minus alpha (this function is responsible for returning the reduced
   // --- cost, which includes alpha).
   // ---
   // --- NOTE, redCost does not include alpha as sent in
   // ---

   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << " <---- solveRelaxed() ---- ";
              (*m_osLog) << "\n";
              );

   m_stats.timerOther1.restart();
  
   //double varRedCost, varOrigCost;
   //bool doPush;
   DecompStat rc;

   //if returns something use it, otherwise don't??

   //if don't have CBC or CPX, this must be defined
  
   m_stats.timerOther2.restart();
   rc = m_app->APPsolveRelaxed(whichModel,
                               redCostX,
                               origCost,
                               alpha,
                               n_origCols,
                               checkRC,
                               checkDup,
                               subprobSI,
                               vars);
   m_stats.thisSolveRelaxApp.push_back(m_stats.timerOther2.timeElapsed());

#if defined (__DECOMP_IP_CBC__) || (__DECOMP_IP_CPX__)
   if(rc == STAT_UNKNOWN){ //sanity check - run this anyway
      //if(1){ //sanity check - run this anyway
      assert(subprobSI);

      // ---
      // --- reset the objective to reduced cost
      // ---
      //TODO: do once, not every time
      int * cinds = new int[n_origCols];
      UtilIotaN(cinds, n_origCols, 0);
      subprobSI->setObjCoeffSet(cinds, cinds + n_origCols, redCostX);
    
      // ---
      // --- solve: min cx, s.t. A'x >= b', x in Z ([A,b] is in modelRelax.M)
      // ---
    
      //TODO: get best N feasible solutions? is this possible with CBC
      //this can help if this is the primary form of col-gen


      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 50,
                 string fileName = "subProbInitCbc";
                 printCurrentProblem(subprobSI, fileName);
                 );
    
#ifdef __DECOMP_IP_CBC__
      //?? does this make a copy? or just pass pointers?
      CbcModel cbc(*subprobSI);
      CbcStrategyDefault cbcStrategyDefault;
      cbc.setStrategy(cbcStrategyDefault);
      cbc.messageHandler()->setLogLevel(m_app->m_param.LogLpLevel);
      //this is NOT a good version of CBC
      cbc.branchAndBound();
      assert(cbc.isProvenOptimal()); //?? NOCOLS cut?
      const double * solution = cbc.getColSolution();
#endif
    
#ifdef __DECOMP_IP_CPX__ 
      subprobSI->branchAndBound(); //this is an instance of LP solver not IP
      const double * solution = subprobSI->getColSolution();
#endif
  
      // THINK: we really don't want to force the user to create vars
      // and check rc and obj, etc... but they might know how to be smart
      // and produce more than one, etc... THINK
    
      //
      // --- create a DecompVar (shat) from the optimal solution
      // --- THINK: assumes 01 everywhere? 
      //
      vector<int>    ind;
      vector<double> els;
      int c;
      int doPush;
      double varRedCost  = 0.0; //stupid - == obj ??
      double varOrigCost = 0.0;
      for(c = 0; c < n_origCols; c++){
         if(fabs(solution[c]) > m_app->m_param.TolZero){
            //printf("\nsol[%d]: %g, redCostX: %g", c, solution[c], redCostX[c]);
        
            //THINK: need to at least be integral?? not if doing DW-LP
        
            ind.push_back(c);
            els.push_back(solution[c]);
        
            //the reduced cost of shat: (c-uA").s
            varRedCost  += redCostX[c] * solution[c];
        
            //the original cost of shat: c.s
            varOrigCost += origCost[c] * solution[c]; 
         }
      }

#ifdef __DECOMP_IP_CBC__
      assert(UtilIsZero(cbc.getObjValue() - varRedCost, m_app->m_param.TolZero));
#else
      assert(UtilIsZero(subprobSI->getObjValue() - varRedCost, m_app->m_param.TolZero));
#endif

      varRedCost += alpha;//RC = c-uA''s - alpha    
      UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                 (*m_osLog) << "\nalpha      = " << alpha;
                 (*m_osLog) << "\nvarRedCost = " << varRedCost;
                 );
    
      // ---
      // --- if checkRC, then only accept variables with negative reduced cost
      // --- if checkDup, check for duplicate variables
      // ---

      //NOTE: empty column is possible - essentially - not really empty
      //once you add the convexity coefficient

      if(ind.size() == 0)
         printf("\nNULL COLUMN");
    
      DecompVar * var = new DecompVar(ind, els, varRedCost, varOrigCost);
    
      doPush = true;
      if(checkRC && varRedCost > -1.e-10) //THINK: dualTol?
         doPush = false;
      else if(checkDup){
         DecompVarList::iterator it;
         for(it = vars.begin(); it != vars.end(); it++){
            if((*it)->isEquivalent(*var)){
               UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
                          (*m_osLog) << "\nDuplicate variable, not adding.";
                          );
               doPush = false;
               break;
            }
         }
      }

      //just use as sanity check doPush = 0
      //doPush = 0;
      if(doPush){
         UTIL_DEBUG(m_app->m_param.LogDebugLevel, 5,
                    var->print();
                    );
         vars.push_back(var);
      }
      else
         UTIL_DELPTR(var);
    
    
      //
      // --- free the local memory
      //
      UTIL_DELARR(cinds);
   }
#endif  
   UTIL_DEBUG(m_app->m_param.LogDebugLevel, 3,
              (*m_osLog) << m_classTag << "  ---- solveRelaxed() ----> ";
              );
   m_stats.thisSolveRelax.push_back(m_stats.timerOther1.timeElapsed()); 
  
   return STAT_UNKNOWN;
}

