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


//thinking about design - PC/C - most everything is driven by OSI methods
//but RC is not, this actually doesn't need OSI at all - either use OsiNull
//or, split up there types...

//DecompAlgo
//   now will have a data member DecompInterface * interface
//   which is allocated to be a DecompOsi or DecompNull depending
//   on which algo we are talking about
// this might also be a way to use Vol directly....


//DecompInterface - base class
//  DecompOsi : public DecompInterface -->
//     wrapper for all OSI methods (with OSI)
//  DecompNull : public DecompInterface -->
//     wrapper for all OSI methods in RC (no OSI)



//===========================================================================//
#ifndef DecompAlgo_h_
#define DecompAlgo_h_

//===========================================================================//
/**
 * \class DecompAlgo
 * \brief Base class for DECOMP algorithms.
 *
 */
//===========================================================================//

//===========================================================================//
#include "Decomp.h"
#include "DecompApp.h"
#include "DecompParam.h"
#include "DecompStats.h"
#include "DecompVarPool.h"
#include "DecompCutPool.h"
#include "DecompMemPool.h"
#include "DecompSolution.h"
#include "DecompAlgoCGL.h"
#include "OsiClpSolverInterface.hpp"
class OsiSolverInterface;
class DecompConstraintSet;
class DecompSolverResult;

//===========================================================================//
class DecompAlgo {

protected:

   //----------------------------------------------------------------------//
   /**
    * @name Data.
    * @{
    */
   //----------------------------------------------------------------------//
   /**
    * Store the name of the class (for logging/debugging) - "who am I?"
    */
   string m_classTag;

   /**
    * Decomp Parameters.
    */
   DecompParam m_param;
   
   /**
    * Type of algorithm for this instance.
    */
   DecompAlgoType m_algo;

   /**
    * The current algorithm status.
    */
   DecompStatus m_status;

   /**
    * The current algorithm phase.
    */
   DecompPhase m_phase;
   DecompPhase m_phaseLast;//just before done
   
   /**
    * Pointer to current active DECOMP application.    
    */
   DecompApp * m_app;

   /**
    * Storage of statistics for run and node.
    */
   DecompStats     m_stats;
   DecompNodeStats m_nodeStats;

   /**
    * Memory pool used to reduce the number of allocations needed.
    */
   DecompMemPool m_memPool;

   /**
    * Stream for log file (default to stdout).
    */
   ostream * m_osLog;

   DecompAlgoCGL * m_cgl;

   /**
    * Pointer (and label) to current active model core/relax.
    */



   
   vector<double>        m_origColLB;
   vector<double>        m_origColUB;

   /**
    * Solver interface(s) for subproblems (P').
    */
   //vector<OsiSolverInterface*> m_subprobSI;
   
   /**
    * Solver interface(s) for master problem (Q'').
    *   CPM: holds model core (and optionally relaxed)
    *        in original space
    *   PC : holds model core in reformulated space
    */
   OsiSolverInterface * m_masterSI;

   /**
    * Solver interface(s) for entire problem (Q'').
    *   CPM: not used (use m_masterSI)
    *   PC : holds model core (and optionally relaxed)
    *        in original space - used for CGL cuts
    */
   OsiClpSolverInterface * m_cutgenSI;
   int                     m_cutgenObjCutInd;
   OsiSolverInterface    * m_auxSI;


   const double                      * m_objective;
   DecompAlgoModel                     m_modelCore;
   map<int, DecompAlgoModel>           m_modelRelax;
   map<int, vector<DecompAlgoModel> >  m_modelRelaxNest;


   /**
    * Containers for variables (current and pool).
    */
   DecompVarList m_vars;
   DecompVarPool m_varpool;

   /**
    * Containers for cuts (current and pool).
    */
   DecompCutList m_cuts;
   DecompCutPool m_cutpool;

   /**
    * Storage for current solution (in x-space).
    */
   double * m_xhat;

   /**
    * User-defined cutoff (global UB) for B&B fathoming and LR.
    *    This does not imply a feasible IP solution, just a bound.
    */
   double   m_cutoffUB;
   //THINK - use solution pool
   vector<DecompSolution*>   m_xhatIPFeas;
   DecompSolution          * m_xhatIPBest;


   //for cpx
   vector<double> colSolution;

   bool m_isColGenExact;

   UtilParameters * m_utilParam;

   int m_numConvexCon;

   //for round robin
   int m_rrLastBlock;
   int m_rrIterSinceAll;

   //
   int m_nArtCols;

   //all these are related to master LP - make object
   int m_nRowsOrig;
   int m_nRowsBranch;
   int m_nRowsConvex;
   int m_nRowsCuts;
   vector<DecompRowType> m_masterRowType;
   vector<DecompColType> m_masterColType;
   vector<int>           m_masterArtCols;

   //to enforce in subproblems
   double * m_colLBNode;
   double * m_colUBNode;

   int      m_compressColsLastPrice;
   int      m_compressColsLastNumCols;

   double         m_relGap;
   DecompAlgoStop m_stopCriteria;
   int            m_colIndexUnique;
   double         m_masterObjLast;//last master obj
   bool           m_objNoChange;


   double       m_stabEpsilon;
   bool         m_useInitLpDuals;
   map<int,int> m_artColIndToRowInd;

   double       m_globalLB;
   double       m_globalUB;

public:   
   /**
    * @}
    */
   //-----------------------------------------------------------------------//
   /**
    * @name Pure virtual functions.
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * Create the master problem (all algorithms must define this function).
    */
   //virtual void createMasterProblem(DecompVarList & initVars) = 0;
   virtual void createMasterProblem(DecompVarList & initVars);
   void loadSIFromModel(OsiSolverInterface           * si,
                        //need next 2 args? ever different?
                        //DecompAlgoModel              & modelCore,
                        //map<int, DecompAlgoModel>    & modelRelax,
                        bool                           doInt = false);
      

   /**
    * Compose solution in x-space from current space.
    *  - PC: this recomposes x from lambda
    *  - C : this just copies over LP solution
    */
   //not pure?
   virtual void recomposeSolution(const double * solution,
                                  double       * rsolution);
   /**
    * @}
    */

   //-----------------------------------------------------------------------//
   /**
    * @name Virtual functions.
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * The main DECOMP process loop for a node.
    */
   virtual DecompStatus processNode(const int nodeIndex   = 0,
                                    const double globalLB = -DecompInf,
                                    const double globalUB =  DecompInf); 
  
   /**
    * Generate initial variables for master problem (PC/DC/RC).
    *   - in CPM, this does nothing
    */
   //THINK: belongs in base? PC or... 
   virtual int generateInitVars(DecompVarList & initVars);

   /**
    * Update of the solution vectors (primal and/or dual).
    */
   virtual DecompStatus 
   solutionUpdate(const DecompPhase phase,
                  const int         maxInnerIter = COIN_INT_MAX,
                  const int         maxOuterIter = COIN_INT_MAX);
   
   /**
    * Update of the phase for process loop.
    */   
   virtual void phaseUpdate(DecompPhase  & phase,
			    DecompStatus & status);
   
   /**
    * Run the initial phase for processing node.
    */      
   virtual void phaseInit(DecompPhase & phase){
      if(getNodeIndex() == 0)
	 phase = PHASE_PRICE1;
   }

   /**
    * Run the done phase for processing node.
    */      
   virtual void phaseDone(){};

   /**
    * Calculate the current LB and update best/history.
    */
   virtual bool updateObjBoundLB(const double mostNegRC = -DecompBigNum);


   virtual void solutionUpdateAsIP(){}

   virtual int adjustColumnsEffCnt(){return DecompStatOk;};
   virtual int compressColumns    (){return DecompStatOk;};
   /**
    * @}
    */
   bool isGapTight(){
      //TODO: make param
      double tightGap = m_param.MasterGapLimit;
      if(m_param.LogDebugLevel >= 2){
	 (*m_osLog) << "DW GAP = " << UtilDblToStr(m_relGap) << "\n";
      }
      if(m_relGap <= tightGap)
	 return true;
      else
	 return false;
   }



   //................................
   //TODO............................
   //................................




   virtual bool isDone() {
      if((m_nodeStats.cutsThisCall + m_nodeStats.varsThisCall) > 0)
	 return false;
      else
	 return true;
   }
   
   //TODO: should move out to PC
   //THINK - helper func?, or specific to PC - right? as is genInit
   vector<double*> getDualRays(int maxNumRays);
   //virtual int generateVarsInf(DecompVarList    & newVars, 
   //		       double           & mostNegReducedCost);
   virtual int generateVarsFea(DecompVarList    & newVars, 
			       double           & mostNegReducedCost);

   virtual int generateVars(const DecompStatus   stat,
			    DecompVarList    & newVars, 
			    double           & mostNegReducedCost);
   virtual int generateCuts(double        * xhat,
			    DecompCutList & newCuts);

   virtual void addVarsToPool(DecompVarList & newVars);
   virtual void addVarsFromPool();
   virtual void addCutsToPool(const double  *  x,
			      DecompCutList & newCuts,
			      int           & m_cutsThisCall);
   virtual int addCutsFromPool();

   bool isIPFeasible(const double * x,
		     const double   feasVarTol = 1.0e-6,  //0.0001%
                     const double   feasConTol = 1.0e-5,  //0.001%
		     const double   intTol     = 1.0e-5); //0.001%

   bool isLPFeasible(const double * x,
		     const double   feasVarTol = 1.0e-6,  //0.001%
                     const double   feasConTol = 1.0e-5); //0.01%

   //fugly
   DecompStatus solveRelaxed(const double        * redCostX,
			     const double        * origCost,
			     const double          alpha,
			     const int             n_origCols,
			     const bool            isNested,
                             DecompAlgoModel     & algoModel,
			     DecompSolverResult  * solveResult,
			     list<DecompVar*>    & vars);
   
   
   inline void appendVars(DecompVar * var){
      m_vars.push_back(var);
   }
   inline void appendVars(DecompVarList & varList){
      copy(varList.begin(), varList.end(), back_inserter(m_vars));
   }
   virtual void setMasterBounds(const double * lbs,
				const double * ubs);

   //here just for RC... ugh...
   //this is where a full working OsiNull might make life easier.
   virtual const double * getRowPrice() const {
      return m_masterSI->getRowPrice();
   }



   int chooseBranchVar(int    & branchedOnIndex,
		       double & branchedOnValue);


   //-----------------------------------------------------------------------//
   /**
    * @name Helper functions.
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * Initial setup of algorithm structures and solver interfaces.
    */
   void initSetup(UtilParameters * utilParam,
		  string         & sectionParam);
   void getModelsFromApp();
   void createOsiSubProblem(DecompAlgoModel & algoModel);
   
   /**
    *
    */
   void coreMatrixAppendColBounds();
   void checkMasterDualObj();
   bool  checkPointFeasible(const DecompConstraintSet * modelCore,
			   const double              * x);
   bool isDualRayInfProof(const double           * dualRay,
			  const CoinPackedMatrix * rowMatrix,
			  const double           * colLB,
			  const double           * colUB,
			  const double           * rowRhs,
			  ostream                * os);
   bool isDualRayInfProofCpx(const double           * dualRay,
			     const CoinPackedMatrix * rowMatrix,
			     const double           * colLB,
			     const double           * colUB,
			     const double           * rowRhs,
			     ostream                * os);

   void printBasisInfo(OsiSolverInterface * si,
		       ostream            * os);


   /**
    *
    */
   void printCurrentProblem(const OsiSolverInterface * si,
                            const string               baseName,
                            const int                  nodeIndex,
                            const int                  cutPass,
                            const int                  pricePass,
                            const int                  blockId    = -1,
                            const bool                 printMps   = false,
                            const bool                 printLp    = true);
   /**
    *
    */
   void printCurrentProblem(const OsiSolverInterface * si,
                            const string               fileName,
                            const bool                 printMps   = true,
                            const bool                 printLp    = true);
   
   /**
    *
    */
   void printVars(ostream * os);
   void printCuts(ostream * os);

   /**
    *
    */
   void createFullMps(const string fileName);

   /**
    *
    */
   virtual void solveDirect(int                  timeLimit = COIN_INT_MAX,
			    DecompSolverResult * result    = NULL){};



   void masterMatrixAddArtCol(CoinPackedMatrix * masterM,
                              char               LorG,
                              int                rowIndex,
                              int                colIndex,
                              DecompColType      colType,
                              double           & colLB,
                              double           & colUB,
                              double           & objCoeff);
   virtual void masterMatrixAddArtCols(CoinPackedMatrix * masterM,
                                       double           * colLB,
                                       double           * colUB,
                                       double           * objCoeff,
                                       vector<string>   & colNames,
                                       int                startRow,
                                       int                endRow,
                                       char               origOrBranch);
   void masterPhaseItoII();
   void masterPhaseIItoI();

   bool isMasterColStructural(const int index) const {
      return (m_masterColType[index] == DecompCol_Structural ||
              m_masterColType[index] == DecompCol_Structural_NoDelete);
   }
   bool isMasterColArtificial(const int index) const {
      return (m_masterColType[index] == DecompCol_ArtForRowL    ||
              m_masterColType[index] == DecompCol_ArtForRowG    ||
              m_masterColType[index] == DecompCol_ArtForBranchL ||
              m_masterColType[index] == DecompCol_ArtForBranchG ||
              m_masterColType[index] == DecompCol_ArtForCutL    ||
              m_masterColType[index] == DecompCol_ArtForCutG);
   }
   


   /**
    * @}
    */
   
   
   //-----------------------------------------------------------------------//
   /**
    * @name Set/get methods.
    * @{
    */
   //-----------------------------------------------------------------------//
   inline const double * getColLBNode() const {
      return m_colLBNode;
   }
   inline const double * getColUBNode() const {
      return m_colUBNode;
   }
   //inline OsiSolverInterface * getSubProbSI(int b){
   // return m_subprobSI[b];
   //}

   inline DecompStats & getStats() {
      return m_stats;
   }

   inline const double * getOrigObjective() const {
     return m_app->m_objective;
   }
   inline const DecompAlgoModel & getModelCore() const {
      return m_modelCore;}

   inline DecompAlgoModel & getModelRelax(const int blockId){
      map<int,DecompAlgoModel>::iterator mit;
      mit = m_modelRelax.find(blockId);
      assert(mit != m_modelRelax.end());      
      return (*mit).second;
   }

   
   /**
    * Get a ptr to the current solution (in x-space).
    */
   inline const double * getXhat() const {
      return m_xhat;
   }

   inline void setCutoffUB(const double thisBound) {
      m_cutoffUB = thisBound;
      setObjBoundUB(thisBound);
   }

   //TODO
   inline const DecompSolution * getXhatIPBest() const {
      return m_xhatIPBest;
   }

   inline const vector<DecompSolution*> & getXhatIPFeas() const {
      return m_xhatIPFeas;
   }
   
   inline const double getCutoffUB() const {
      return m_cutoffUB;
   }
   
   inline DecompStats & getDecompStats(){
      return m_stats;
   }

   inline const DecompParam & getDecompParam() const{
      return m_param;
   }

   inline const DecompApp * getDecompApp() const{
      return m_app;
   }
   inline DecompApp * getDecompAppMutable() {
      return m_app;
   }

   inline const int getNodeIndex() const {
      return m_nodeStats.nodeIndex;
   }

   inline const int getCutCallsTotal() const {
      return m_nodeStats.cutCallsTotal;
   }

   inline const int getPriceCallsTotal() const {
      return m_nodeStats.priceCallsTotal;
   }

#if 1
   inline const double * getMasterColSolution() const {
      return &colSolution[0];
   }

   inline double getMasterObjValue() const {
      //don't believe this, could have just added a col
      //int nc = m_masterSI->getNumCols();

      //be careful that this is always using the PhaseII obj
      
      //double * objCoef   = m_app->m_model.objCoeff; //THINK

      int nc = static_cast<int>(colSolution.size());
      const double * objCoef = m_masterSI->getObjCoefficients();
      const double * colSol  = getMasterColSolution();  
      double retVal = 0.0;
      for ( int i=0 ; i<nc ; i++ ){
	 retVal += objCoef[i]*colSol[i];
      }
      return retVal;
   }

#else
   inline const double * getMasterColSolution() const {
      return m_masterSI->getColSolution();
   }
   inline double getMasterObjValue() const {
      return m_masterSI->getObjValue(); //if we use base version
   }
#endif


   inline const int getStopCriteria() const {
      return m_stopCriteria;
   }

   /** 
    * Get the current best LB.
    */
   inline const double getObjBestBoundLB() const {
      return m_nodeStats.objBest.first;
   }

   /**
    * Get the current best UB.
    */
   inline const double getObjBestBoundUB() const {
      return m_nodeStats.objBest.second;
   }

   /**
    * Set the current LB and update best/history.
    */
   inline void setObjBoundLB(const double thisBound){
      if(thisBound > m_nodeStats.objBest.first){
	 m_nodeStats.objBest.first = thisBound;
      }

      DecompObjBound objBound;
      objBound.lbOrUb    = 0;
      objBound.cutPass   = m_nodeStats.cutCallsTotal;
      objBound.pricePass = m_nodeStats.priceCallsTotal;      
      objBound.thisBound = thisBound;
      objBound.bestBound = m_nodeStats.objBest.first;
#ifdef UTIL_USE_TIMERS
      objBound.timeStamp = globalTimer.getRealTime();
#else
      objBound.timeStamp = -1;
#endif
      
      m_nodeStats.objHistoryLB.push_back(objBound);
   }
   
   /**
    * Set the current UB and update best/history.
    */
   inline void setObjBoundUB(const double thisBound){
      if(thisBound < m_nodeStats.objBest.second){
         //printf("thisBound= %g objBest= %g\n",
         //     thisBound, m_nodeStats.objBest.second);
	 UTIL_MSG(m_app->m_param.LogDebugLevel, 3,
		  (*m_osLog) << "New Global UB = " 
		  << UtilDblToStr(thisBound) << endl;);
	 m_nodeStats.objBest.second = thisBound;
         
      }

      DecompObjBound objBound;
      objBound.lbOrUb    = 1;
      objBound.cutPass   = m_nodeStats.cutCallsTotal;
      objBound.pricePass = m_nodeStats.priceCallsTotal;      
      objBound.thisBound = thisBound;
      objBound.bestBound = m_nodeStats.objBest.second;
#ifdef UTIL_USE_TIMERS
      objBound.timeStamp = globalTimer.getRealTime();
#else
      objBound.timeStamp = -1;
#endif      
      m_nodeStats.objHistoryUB.push_back(objBound);
   }




   bool isTailoffLB(const int    changeLen      = 10,
		    const double changePerLimit = 0.1);


   inline int getNumRowType(DecompRowType rowType){
      int   nRowsType = 0;
      vector<DecompRowType>::iterator vi;
      for(vi = m_masterRowType.begin(); vi != m_masterRowType.end(); vi++){
	 if(*vi == rowType)
	    nRowsType++;
      }
      return nRowsType;
   }
   
   void checkBlocksColumns();


   /**
    * @}
    */

   //TODO:
   //be careful here that we don't stop due to mLB>=m_UB in the case where
   //user gives optimal UB as cutoff, but we don't yet have integral solution
   


   
   //-----------------------------------------------------------------------//
   /**
    * @name Constructors and destructor.
    * @{
    */
   //-----------------------------------------------------------------------//
public:
   /**
    * Default constructors.
    */   
   DecompAlgo(const DecompAlgoType   algo,
              DecompApp            * app,
              UtilParameters       * utilParam):
      m_classTag   ("D-ALGO"),
      m_param      (),
      m_algo       (algo),
      m_status     (STAT_UNKNOWN),
      m_phase      (PHASE_UNKNOWN),
      m_phaseLast  (PHASE_UNKNOWN),
      m_app        (app),
      m_stats      (),
      m_nodeStats  (),
      m_memPool    (),
      m_osLog      (&cout),
      m_cgl          (0),
      m_origColLB  (),
      m_origColUB  (),
      m_masterSI   (0),
      m_cutgenSI   (NULL),
      m_cutgenObjCutInd(-1),
      m_auxSI      (NULL),
      m_vars       (),
      m_varpool    (),
      m_cuts       (),
      m_cutpool    (),
      m_xhat       (0),
      m_cutoffUB   (DecompInf),
      m_xhatIPFeas (),
      m_xhatIPBest (NULL),
      m_isColGenExact(false),
      m_utilParam    (utilParam),
      m_numConvexCon (1),
      m_rrLastBlock (-1),
      m_rrIterSinceAll(0),

      m_colLBNode(NULL),
      m_colUBNode(NULL),
      m_relGap(DecompInf),
      m_stopCriteria(DecompStopNo),
      m_masterObjLast(DecompInf)
   {
      m_app->m_decompAlgo = this;
   }
   

   /**
    * Destructor.
    */
   virtual ~DecompAlgo(){
      //UtilDeleteVectorPtr(m_subprobSI);
      UTIL_DELPTR(m_masterSI);
      UTIL_DELPTR(m_cutgenSI);
      UTIL_DELPTR(m_auxSI);
      UTIL_DELARR(m_xhat);
      UTIL_DELPTR(m_cgl);
      UtilDeleteVectorPtr(m_xhatIPFeas);
      UtilDeleteListPtr(m_vars);
      UtilDeleteListPtr(m_cuts);
      UTIL_DELARR(m_colLBNode);
      UTIL_DELARR(m_colUBNode);
   }
   /**
    * @}
    */



};

#endif
