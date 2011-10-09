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

#ifndef DECOMP_ALGO_INCLUDED
#define DECOMP_ALGO_INCLUDED

#define EXPLICIT_BOUNDS

#include "DecompConstraintSet.h"
#include "DecompSolution.h"
#include "DecompMemPool.h"
#include "DecompVarPool.h"
#include "DecompCutPool.h"
#include "DecompStats.h"
//class DecompStats;
class DecompApp;
class ClpModel;


#include "OsiSolverInterface.hpp"



// --------------------------------------------------------------------- //
class DecompAlgo{
 private:
   DecompAlgo(const DecompAlgo &);
   DecompAlgo & operator=(const DecompAlgo &);
  
 private:
   static const char * m_classTag;
   static const char * versionTag;
   DecompAlgoType      m_algo;     //"who am i"
  
 protected:
   DecompMemPool         m_auxMemPool;
   DecompApp           * m_app;

   int                   m_nodeIndex;
   int                   m_whichModel;
   int                   m_whichCoreModel;
   int                   m_priceCallsRound;
   int                   m_priceCallsTotal;
   int                   m_cutCallsRound;
   int                   m_cutCallsTotal;
   int                   m_varsThisRound;
   int                   m_cutsThisRound;
   int                   m_varsThisCall;
   int                   m_cutsThisCall;

   int                   m_isTightenAlgo;

   ostream            *  m_osLog;
  
   //THINK: naming outer, inner?? - will work across all?

   //solver interface for subproblem (P')  
   map<int, OsiSolverInterface*>   m_subprobSI;

   //solver interface for master (Q")
   OsiSolverInterface          * m_masterSI;

#ifdef __DECOMP_LP_CLP__
   //THINK: so we don't freecached when getModelPtr - do once... 
   //makes it so that we are stuck with CLP - also need to support CPX
   ClpModel                    * m_masterCLP;
#endif


   //m_vars and m_cuts have no destructor, but they are list of ptrs
   //the memory they point to needs to be deleted.
  
   //member of app or member of algo?
   DecompVarList           m_vars;          //list of vars added to master
   DecompVarList           m_initVars;
   DecompCutList           m_cuts;
   DecompVarPool           m_varpool;
   DecompCutPool           m_cutpool;

   double                  m_tlb;
   double                  m_tub;
   double                  m_bestUpperBound;//known opt

   double                  * m_xhat;
   
   //think about do we want a solution pool?
   //for now, just keep the best and toss the rest
   vector<DecompSolution*>   m_xhatIPFeas;
   DecompSolution          * m_xhatIPBest;

   int m_numOrigCols; //DECOMP

   vector< vector<double> > m_optPoint;

   //stabilization techniques
   double * m_piEstimate;
   vector<bool> isStab;

 public:
   //these are just pointers into app
   DecompConstraintSet * m_modelCore;
   DecompConstraintSet * m_modelRelax;

   DecompStats           m_stats;

 public:
   //helper functions
   OsiSolverInterface * initSolverInterface(); 
  
   //TODO: functions in alphabetical order??
   void   startupLog();
   void   initSetup(int whichModel = 1);
  
   void   tighten(int whichModel);
   inline double   getTrueLowerBound(){
      return m_tlb;
   }
   inline double   getTrueUpperBound(){
      return m_tub;
   }
   
   //TODO: should not call these LB UB
   void   setTrueLowerBound(const double mostNegReducedCost);
   void   setTrueUpperBound(const double ub){
      m_tub = ub;
   };
   double calcConstant(const int      m, 
                       const double * u);

   //helper functions - DecompDebug.cpp
   bool isDualRayInfProof(const double           * dualRay,
                          const CoinPackedMatrix * rowMatrix,
                          const double           * colLB,
                          const double           * colUB,
                          const double           * rowRhs,
                          ostream                * os     = 0);
   void printBasisInfo(OsiSolverInterface * si,
                       ostream            * os);
   void printCurrentProblem(const OsiSolverInterface * si,
                            const string               baseName,
                            const int                  nodeIndex,
                            const int                  cutPass,
                            const int                  pricePass,
                            const bool                 printMps   = true,
                            const bool                 printLp    = true);
   void printCurrentProblem(const OsiSolverInterface * si,
                            const string               fileName,
                            const bool                 printMps   = true,
                            const bool                 printLp    = true);

   void   printVars(ostream * os = &cout);
   void   printCuts(ostream * os = &cout);

   void   solveBruteForce();
   void   createFullMps(const string filename);
   vector<double*> getDualRays(int maxNumRays);

   inline void setApp(DecompApp * app){
      m_app = app;
   }
  
  
   inline void setBestUpperBound(const double bestUpperBound){
      m_bestUpperBound = bestUpperBound;
   }

   DecompStat solveRelaxed(const int             whichModel,
                           const double        * redCostX,
                           const double        * origCost,
                           const double          alpha,
                           const int             n_origCols,
                           const bool            checkRC,
                           const bool            checkDup,
                           OsiSolverInterface  * m_subprobSI,
                           list<DecompVar*>    & vars);

 public:
   //pure virtual methods




   virtual void createMasterProblem(DecompVarList & initVars) = 0;

   //possible base is enough here too... think at least for PC and C 
   //seem to be ok

   //just base is enough for PC and C here
   virtual DecompStat solutionUpdate(const DecompPhase phase,
                                     const int         maxInnerIter,
                                     const int         maxOuterIter);

   //just base is enough?
   virtual DecompPhase phaseUpdate(const DecompPhase   phase,
                                   const DecompStat    stat);
  
 public:
   //virtual methods
   //all this mess, just so that we don't use m_masterSI for RC??
   //doesn't make alot of sense to me... THINK
  
   //just feed u vector in rowPrice in OSI???

   virtual void setMasterBounds(const double * lbs,
                                const double * ubs);

   OsiSolverInterface * getMasterSolverInterface(){
      return m_masterSI;
   }
  
   virtual const double * getRowPrice() const {
      return m_masterSI->getRowPrice();
   }
  
   inline const double * getX(){
      return m_xhat;
   }

   inline DecompApp * getApp(){
      return m_app;
   }

   inline const DecompSolution * getXhatIPBest(){
      return m_xhatIPBest;
   }

#if 1
   virtual const double * getRightHandSide() const {
      return &m_modelCore->rowRhs[0];//THINK
      //return m_masterSI->getRightHandSide();
   }
   virtual const char * getRowSense() const {
      return &m_modelCore->rowSense[0];//THINK - should be SI?
      //return m_masterSI->getRowSense();
   }
#endif

   int heuristics(const double            * xhat,
                  vector<DecompSolution*> & xhatIPFeas);

   virtual int generateVars(const DecompStat   stat,
                            DecompVarList    & newVars, 
                            double           & mostNegReducedCost);
   virtual int generateCuts(DecompCutList & newCuts);

   //different name? set x vector? and maybe base just memcpy
   //while PC does recompose then copies in? 
   virtual void recomposeSolution(const double * solution,
                                  double       * rsolution){
      printf("\nbase recomposeSolution does nothing.");
   }; 

   virtual int  generateInitVars(DecompVarList & initVars);
   virtual bool isDone() { return false; };
  

   //will never get called in C, yet PC specific, do nothing in base
   //and move to PC?
   void addVarsToPool(DecompVarList & newVars);
   void addVarsFromPool();  


   bool isIPFeasible(const double * x,
                     const double   feasTol = 1.0e-4,
                     const double   intTol  = 1.0e-4);
   bool isLPFeasible(const double * x,
                     const double   feasTol = 1.0e-4);
  
  
   //different for C and PC, pure virtual for now
   virtual void addCutsToPool(const double  *  x,
                              DecompCutList & newCuts,
                              int           & n_newCuts);
   virtual int addCutsFromPool();

   //DecompBranch.cpp
   int chooseBranchVar(int    & branchedOnIndex,
                       double & branchedOnValue);
   virtual int branch(int    branchedOnIndex,
                      double branchedOnValue);
   DecompStat processNode(const int nodeIndex   = 0,
                          const double globalLB = -DecompInf,
                          const double globalUB =  DecompInf); 
   //process node? or just process?


   /*helper functions*/
   inline void appendVars(DecompVar * var){
      m_vars.push_back(var);
   }
   inline void appendVars(DecompVarList & varList){
      copy(varList.begin(), varList.end(), back_inserter(m_vars));
#if 0
      if(m_vars.empty()){
         m_vars = varList;
      }
      else
         m_vars.insert(m_vars.end(), varList.begin(), varList.end());
#endif
   }

   
  
 public:
   DecompAlgo(const DecompAlgoType   algo,
              DecompApp            * app) :
      m_algo(algo),
      m_auxMemPool(),
      m_app(app),
      m_whichModel(-1),
      m_whichCoreModel(-1),
      m_nodeIndex(0),
      m_priceCallsRound(0),
      m_priceCallsTotal(0),
      m_cutCallsRound(0),
      m_cutCallsTotal(0),
      m_varsThisRound(0),
      m_cutsThisRound(0),
      m_varsThisCall(0),
      m_cutsThisCall(0),
      m_isTightenAlgo(0),
      m_osLog(&cout),//TODO
      m_masterSI(0),
      m_vars(),
      m_initVars(),
      m_cuts(),
      m_varpool(),
      m_cutpool(),
      m_tlb(-DecompInf),
      m_tub( DecompInf),
      m_bestUpperBound(DecompInf),
      m_xhat(0),
      m_xhatIPBest(0),
      m_piEstimate(0) 
      {
         (*m_osLog).setf(ios::fixed);
      };
  
   virtual ~DecompAlgo(){
      map<int, OsiSolverInterface*>::iterator it;
      for(it = m_subprobSI.begin(); 
          it != m_subprobSI.end(); it++){
         if(it->second)
            UTIL_DELPTR(it->second);      
      }
      UTIL_DELPTR(m_masterSI);
      UTIL_DELARR(m_xhat);
      UtilDeleteVectorPtr(m_xhatIPFeas);
    
      //THINK - if was class, the std::list<T*> would be cleaner?
      //let T's destructor do the work? the way cutpool manages that?
      //should m_cuts and m_vars be pools anyway? vector to waiting row
      //or list of waiting row pts?
      UtilDeleteListPtr(m_cuts);
      UtilDeleteListPtr(m_vars);

      //THINK - m_initVars always (?) gets pushed into m_vars -
      //so if we free here we'll double free?
      //UtilDeleteListPtr(m_initVars);



      //?? who deletes m_var, m_cuts????
   };
};

#endif
        
