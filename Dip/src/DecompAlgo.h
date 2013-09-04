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

//#define DECOMP_MASTERONLY_DIRECT

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
#include "AlpsDecompTreeNode.h"
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
   std::string m_classTag;

   /**
   // DIP is distributed under the Eclipse Public License as part of the        //
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
   DecompPhase m_phaseForce;

   /**
    * Pointer to current active DECOMP application.
    */
   DecompApp* m_app;

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
   std::ostream* m_osLog;

   DecompAlgoCGL* m_cgl;

   /**
    * Pointer (and label) to current active model core/relax.
    */




   std::vector<double>        m_origColLB;
   std::vector<double>        m_origColUB;

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
   OsiSolverInterface* m_masterSI;

   /**
    * Solver interface(s) for entire problem (Q'').
    *   CPM: not used (use m_masterSI)
    *   PC : holds model core (and optionally relaxed)
    *        in original space - used for CGL cuts
    */
   OsiClpSolverInterface* m_cutgenSI;
   int                     m_cutgenObjCutInd;
   OsiSolverInterface*     m_auxSI;


   const double*                       m_objective;
   DecompAlgoModel                     m_modelCore;
   std::map<int, DecompAlgoModel>           m_modelRelax;
   std::map<int, std::vector<DecompAlgoModel> >  m_modelRelaxNest;


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
   double* m_xhat;

   /**
    * User-defined cutoff (global UB) for B&B fathoming and LR.
    *    This does not imply a feasible IP solution, just a bound.
    */
   double   m_cutoffUB;
   //THINK - use solution pool
   std::vector<DecompSolution*>   m_xhatIPFeas;
   DecompSolution*           m_xhatIPBest;


   //for cpx
   std::vector<double> m_primSolution;
   std::vector<double> m_dualSolution;

   bool m_isColGenExact;

   UtilParameters* m_utilParam;

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
   std::vector<DecompRowType> m_masterRowType;
   std::vector<DecompColType> m_masterColType;
   std::vector<int>           m_masterArtCols;

   //to enforce in subproblems
   double* m_colLBNode;
   double* m_colUBNode;

   int      m_compressColsLastPrice;
   int      m_compressColsLastNumCols;

   /**
    * Current node gap (bestUB-bestLB)/bestLB.
    */
   double         m_relGap;

   DecompAlgoStop m_stopCriteria;
   int            m_colIndexUnique;
   double         m_masterObjLast;//last master obj
   bool           m_objNoChange;


   double       m_stabEpsilon;
   bool         m_useInitLpDuals;
   std::map<int, int> m_artColIndToRowInd;

   double       m_globalLB;
   double       m_globalUB;

   std::vector<double> m_phaseIObj;

   int          m_function;//calling function
   bool         m_firstPhase2Call;
   bool         m_isStrongBranch;

   const AlpsDecompTreeNode* m_curNode;

#ifdef DECOMP_MASTERONLY_DIRECT
   //NOTE:
   // this should be found by framework
   //   for first pass, have it set by user (MILPBlock)
   vector<int>  m_masterOnlyCols;
   //vector<bool> m_isColMasterOnly;
   /**
    *  Map from original index to master index for master-only vars.
    */
   map<int, int> m_masterOnlyColsMap;
#endif

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
   virtual void createMasterProblem(DecompVarList& initVars);
   void loadSIFromModel(OsiSolverInterface*            si,
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
   virtual void recomposeSolution(const double* solution,
                                  double*        rsolution);
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
   virtual DecompStatus processNode(const AlpsDecompTreeNode* node,
                                    const double globalLB = -DecompInf,
                                    const double globalUB =  DecompInf);

   /**
    * Provide the current node the algorithm is solving.
    */

   const AlpsDecompTreeNode* getCurrentNode() const {
      return m_curNode;
   };

   /**
    * Do some information sending after the current node has been processed.
    * Does nothing by default.
    */

   virtual void postProcessNode(DecompStatus decompStatus) {};

   /**
    * Do some information sending after the current node has been branched.
    * Does nothing by default.
    */

   virtual void postProcessBranch(DecompStatus decompStatus) {};

   /**
    * Generate initial variables for master problem (PC/DC/RC).
    *   - in CPM, this does nothing
    */
   //THINK: belongs in base? PC or...
   virtual int generateInitVars(DecompVarList& initVars);

   /**
    * Update of the solution vectors (primal and/or dual).
    */
   virtual DecompStatus
   solutionUpdate(const DecompPhase phase,
                  const bool        resolve = true,
                  const int         maxInnerIter = COIN_INT_MAX,
                  const int         maxOuterIter = COIN_INT_MAX);

   /**
    * Update of the phase for process loop.
    */
   virtual void phaseUpdate(DecompPhase&   phase,
                            DecompStatus& status);

   /**
    * Run the initial phase for processing node.
    */
   virtual void phaseInit(DecompPhase& phase) {
      if (getNodeIndex() == 0) {
         phase = PHASE_PRICE1;
      }
   }

   /**
    * Run the done phase for processing node.
    */
   virtual void phaseDone() {};

   /**
    * Calculate the current LB and update best/history.
    */
   virtual bool updateObjBound(const double mostNegRC = -DecompBigNum);


   virtual void solutionUpdateAsIP() {}

   virtual int adjustColumnsEffCnt() {
      return DecompStatOk;
   };
   virtual int compressColumns    () {
      return DecompStatOk;
   };
   /**
    * @}
    */
   bool isGapTight() {
      //TODO: make param
      double tightGap = m_param.MasterGapLimit;

      //printf("isGapTight m_relGap = %g\n", m_relGap);
      if (m_param.LogDebugLevel >= 2) {
         (*m_osLog) << "DW GAP = " << UtilDblToStr(m_relGap)
                    << " isTight = " << (m_relGap <= tightGap)
                    << "\n";
      }

      if (m_relGap <= tightGap) {
         return true;
      } else {
         return false;
      }
   }



   //................................
   //TODO............................
   //................................




   virtual bool isDone() {
      if ((m_nodeStats.cutsThisCall + m_nodeStats.varsThisCall) > 0) {
         return false;
      } else {
         return true;
      }
   }

   //TODO: should move out to PC
   //THINK - helper func?, or specific to PC - right? as is genInit
   std::vector<double*> getDualRays(int maxNumRays);
   virtual int generateVarsFea(DecompVarList&     newVars,
                               double&            mostNegReducedCost);

   virtual int generateVars(const DecompStatus   stat,
                            DecompVarList&     newVars,
                            double&            mostNegReducedCost);
   virtual int generateCuts(double*         xhat,
                            DecompCutList& newCuts);

   virtual void addVarsToPool(DecompVarList& newVars);
   virtual void addVarsFromPool();
   virtual void addCutsToPool(const double*    x,
                              DecompCutList& newCuts,
                              int&            m_cutsThisCall);
   virtual int addCutsFromPool();

   bool isIPFeasible(const double* x,
                     const bool     isXSparse  = false,
                     const double   feasVarTol = 1.0e-6,  //0.0001%
                     const double   feasConTol = 1.0e-5,  //0.001%
                     const double   intTol     = 1.0e-5); //0.001%

   bool isLPFeasible(const double* x,
                     const bool     isXSparse  = false,
                     const double   feasVarTol = 1.0e-6,  //0.001%
                     const double   feasConTol = 1.0e-5); //0.01%

   //fugly
   DecompStatus solveRelaxed(const double*         redCostX,
                             const double*         origCost,
                             const double          alpha,
                             const int             n_origCols,
                             const bool            isNested,
                             DecompAlgoModel&      algoModel,
                             DecompSolverResult*   solveResult,
                             std::list<DecompVar*>&     vars);


   inline void appendVars(DecompVar* var) {
      m_vars.push_back(var);
   }
   inline void appendVars(DecompVarList& varList) {
      copy(varList.begin(), varList.end(), back_inserter(m_vars));
   }
   virtual void setMasterBounds(const double* lbs,
                                const double* ubs);
   virtual void setSubProbBounds(const double* lbs,
                                 const double* ubs);

   //int chooseBranchVar(int    & branchedOnIndex,
   //	       double & branchedOnValue);
   virtual bool
   chooseBranchSet(std::vector< std::pair<int, double> >& downBranchLb,
                   std::vector< std::pair<int, double> >& downBranchUb,
                   std::vector< std::pair<int, double> >& upBranchLb,
                   std::vector< std::pair<int, double> >& upBranchUb);




   //-----------------------------------------------------------------------//
   /**
    * @name Helper functions.
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * Initial setup of algorithm structures and solver interfaces.
    */
   void initSetup(UtilParameters* utilParam,
                  std::string&          sectionParam);
   void getModelsFromApp();
   void createOsiSubProblem(DecompAlgoModel& algoModel);

   /**
    * Calculate gap: |(ub-lb)|/|lb|
    */
   /*double calculateGap(double boundLB,
   	       double boundUB) const {
      double gap = DecompInf;
      if(boundLB > -DecompInf && boundUB < DecompInf){
    if(boundLB != 0.0)
       gap = fabs(boundUB-boundLB)/fabs(boundLB);
    else
       gap = fabs(boundUB);
      }
      return gap;
      }*/


   /**
    *
    */
   void coreMatrixAppendColBounds();
   void checkMasterDualObj();
   bool  checkPointFeasible(const DecompConstraintSet* modelCore,
                            const double*               x);
   bool isDualRayInfProof(const double*            dualRay,
                          const CoinPackedMatrix* rowMatrix,
                          const double*            colLB,
                          const double*            colUB,
                          const double*            rowRhs,
                          std::ostream*                 os);
   bool isDualRayInfProofCpx(const double*            dualRay,
                             const CoinPackedMatrix* rowMatrix,
                             const double*            colLB,
                             const double*            colUB,
                             const double*            rowRhs,
                             std::ostream*                 os);

   void printBasisInfo(OsiSolverInterface* si,
                       std::ostream*             os);


   /**
    *
    */
   void printCurrentProblemDual(OsiSolverInterface* si,
                                const std::string         baseName,
                                const int            nodeIndex,
                                const int            cutPass,
                                const int            pricePass);

   void printCurrentProblem(const OsiSolverInterface* si,
                            const std::string               baseName,
                            const int                  nodeIndex,
                            const int                  cutPass,
                            const int                  pricePass,
                            const int                  blockId    = -1,
                            const bool                 printMps   = true,//false,
                            const bool                 printLp    = true);
   /**
    *
    */
   void printCurrentProblem(const OsiSolverInterface* si,
                            const std::string               fileName,
                            const bool                 printMps   = true,
                            const bool                 printLp    = true);

   /**
    *
    */
   void printVars(std::ostream* os);
   void printCuts(std::ostream* os);

   /**
    *
    */
   void createFullMps(const std::string fileName);

   /**
    *
    */
   virtual DecompSolverResult*
   solveDirect(const DecompSolution* startSol  = NULL) {
      return NULL;
   }

#ifdef DECOMP_MASTERONLY_DIRECT
   void masterMatrixAddMOCols(CoinPackedMatrix* masterM,
                              double*            colLB,
                              double*            colUB,
                              double*            objCoeff,
                              std::vector<string>&    colNames);
#endif

   void masterMatrixAddArtCol(std::vector<CoinBigIndex>& colBeg,
                              std::vector<int         >& colInd,
                              std::vector<double      >& colVal,
                              char                   LorG,
                              int                    rowIndex,
                              int                    colIndex,
                              DecompColType          colType,
                              double&                colLB,
                              double&                colUB,
                              double&                objCoeff);

   virtual void masterMatrixAddArtCols(CoinPackedMatrix* masterM,
                                       double*            colLB,
                                       double*            colUB,
                                       double*            objCoeff,
                                       std::vector<std::string>&    colNames,
                                       int                startRow,
                                       int                endRow,
                                       DecompRowType      rowType);
   void masterPhaseItoII();
   void masterPhaseIItoI();

   bool isMasterColMasterOnly(const int index) const {
      return (m_masterColType[index] == DecompCol_MasterOnly);
   }
   bool isMasterColStructural(const int index) const {
      return (m_masterColType[index] == DecompCol_Structural ||
              m_masterColType[index] == DecompCol_Structural_NoDelete);
   }
   bool isMasterColArtificial(const int index) const {
      return (m_masterColType[index] == DecompCol_ArtForRowL    ||
              m_masterColType[index] == DecompCol_ArtForRowG    ||
              m_masterColType[index] == DecompCol_ArtForBranchL ||
              m_masterColType[index] == DecompCol_ArtForBranchG ||
              m_masterColType[index] == DecompCol_ArtForConvexL ||
              m_masterColType[index] == DecompCol_ArtForConvexG ||
              m_masterColType[index] == DecompCol_ArtForCutL    ||
              m_masterColType[index] == DecompCol_ArtForCutG);
   }

   void breakOutPartial(const double*   xHat,
                        DecompVarList& newVars,
                        const double    intTol = 1.0e-5);

   /**
    * Create an adjusted dual vector with the duals from the
    * convexity constraints removed.
    */
   void generateVarsAdjustDuals(const double* uOld,
                                double*        uNew);
   /**
    * Calculated reduced cost vector (over vars in compact space)
    * for a given dual vector.
    */
   void generateVarsCalcRedCost(const double* u,
                                double*        redCostX);




   /**
    * @}
    */


   //-----------------------------------------------------------------------//
   /**
    * @name Set/get methods.
    * @{
    */
   //-----------------------------------------------------------------------//
   inline const double* getColLBNode() const {
      return m_colLBNode;
   }
   inline const double* getColUBNode() const {
      return m_colUBNode;
   }
   //inline OsiSolverInterface * getSubProbSI(int b){
   // return m_subprobSI[b];
   //}

   inline DecompStats& getStats() {
      return m_stats;
   }

   inline const double* getOrigObjective() const {
      return m_app->m_objective;
   }
   inline const DecompAlgoModel& getModelCore() const {
      return m_modelCore;
   }

   inline const int getAlgo() const {
      return m_algo;
   }

   inline const DecompParam& getParam() const {
      return m_param;
   }

   inline DecompParam& getMutableParam() {
      return m_param;
   }

   inline OsiSolverInterface* getMasterOSI() {
      return m_masterSI;
   }

   inline DecompAlgoModel& getModelRelax(const int blockId) {
      std::map<int, DecompAlgoModel>::iterator mit;
      mit = m_modelRelax.find(blockId);
      assert(mit != m_modelRelax.end());
      return (*mit).second;
   }


   /**
    * Get a ptr to the current solution (in x-space).
    */
   inline const double* getXhat() const {
      return m_xhat;
   }

   inline void setCutoffUB(const double thisBound) {
      m_cutoffUB = thisBound;
      setObjBoundIP(thisBound);
   }

   //TODO
   inline const DecompSolution* getXhatIPBest() const {
      return m_xhatIPBest;
   }

   inline const std::vector<DecompSolution*>& getXhatIPFeas() const {
      return m_xhatIPFeas;
   }

   inline const double getCutoffUB() const {
      return m_cutoffUB;
   }

   inline DecompStats& getDecompStats() {
      return m_stats;
   }

   inline const DecompParam& getDecompParam() const {
      return m_param;
   }

   inline const DecompApp* getDecompApp() const {
      return m_app;
   }
   inline DecompApp* getDecompAppMutable() {
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

   /**
    * Get current primal solution for master problem.
    */
   inline const double* getMasterPrimalSolution() const {
      return &m_primSolution[0];
   }

   /**
    * Get current dual solution for master problem.
    */
   virtual const double* getMasterDualSolution() const {
      return &m_dualSolution[0];
   }

   /**
    * Adjust the current dual solution for master problem.
    */
   virtual void adjustMasterDualSolution() {};


   inline double getMasterObjValue() const {
      if (!m_masterSI) {
         return -DecompInf;
      }

      //NOTE: be careful that this is always using the PhaseII obj
      int nc = static_cast<int>(m_primSolution.size());
      const double* objCoef = m_masterSI->getObjCoefficients();
      const double* primSol  = getMasterPrimalSolution();
      double retVal = 0.0;

      for ( int i = 0 ; i < nc ; i++ ) {
         retVal += objCoef[i] * primSol[i];
      }

      return retVal;
   }

   inline const int getStopCriteria() const {
      return m_stopCriteria;
   }

   /**
    * Get the current global (integrality) gap.
    */
   inline const double getGlobalGap() const {
      return UtilCalculateGap(m_globalLB, m_globalUB);
   }

   /**
    * Get the current node (integrality) gap.
    */
   inline const double getNodeIPGap() const {
      return UtilCalculateGap(getObjBestBoundLB(), getObjBestBoundUB());
   }

   /**
    * Get the current node (continuous) gap.
    */
   inline const double getNodeLPGap() const {
      int nHistorySize
         = static_cast<int>(m_nodeStats.objHistoryBound.size());

      if (nHistorySize > 0) {
         const DecompObjBound& objBound
            = m_nodeStats.objHistoryBound[nHistorySize - 1];
         return UtilCalculateGap(getObjBestBoundLB(), objBound.thisBoundUB);
      } else {
         return DecompInf;
      }
   }

   /**
    * Get the current best LB.
    */
   inline const double getObjBestBoundLB() const {
      return m_nodeStats.objBest.first;
   }

   /**
    * Set the object to be in strong branching mode.
    */
   inline const void setStrongBranchIter(bool isStrongBranch = true) {
      m_isStrongBranch = isStrongBranch;
   }

   /**
    * Get the current best UB.
    */
   inline const double getObjBestBoundUB() const {
      return m_nodeStats.objBest.second;
   }

   /**
    * Get a specific row type.
    */
   inline const double getMasterRowType(int row) const {
      return m_masterRowType[row];
   }

   /**
    * Set the current continuous bounds and update best/history.
    */
   virtual void setObjBound(const double thisBound,
                            const double thisBoundUB) {
      UtilPrintFuncBegin(m_osLog, m_classTag,
                         "setObjBound()", m_param.LogDebugLevel, 2);

      if (thisBound > m_nodeStats.objBest.first) {
         m_nodeStats.objBest.first = thisBound;

         if (getNodeIndex() == 0) {
            m_globalLB = thisBound;
         }
      }

      DecompObjBound objBound;
      objBound.phase         = m_phase == PHASE_PRICE1 ? 1 : 2;
      objBound.cutPass       = m_nodeStats.cutCallsTotal;
      objBound.pricePass     = m_nodeStats.priceCallsTotal;
      objBound.thisBound     = thisBound;
      objBound.thisBoundUB   = thisBoundUB;
      objBound.bestBound     = m_nodeStats.objBest.first;
      objBound.bestBoundIP   = m_nodeStats.objBest.second;
#ifdef UTIL_USE_TIMERS
      objBound.timeStamp     = globalTimer.getRealTime();
#else
      objBound.timeStamp     = -1;
#endif
      m_nodeStats.objHistoryBound.push_back(objBound);
      UtilPrintFuncEnd(m_osLog, m_classTag,
                       "setObjBound()", m_param.LogDebugLevel, 2);
   }

   /**
    * Set the current integer bound and update best/history.
    */
   virtual inline void setObjBoundIP(const double thisBound) {
      UtilPrintFuncBegin(m_osLog, m_classTag,
                         "setObjBoundIP()", m_param.LogDebugLevel, 2);

      if (thisBound < m_nodeStats.objBest.second) {
         UTIL_MSG(m_app->m_param.LogDebugLevel, 3,
                  (*m_osLog) << "New Global UB = "
                  << UtilDblToStr(thisBound) << std::endl;);
         m_nodeStats.objBest.second = thisBound;
      }

      //---
      //--- copy the last continuous history, adjust the time
      //---
      DecompObjBound   objBoundIP;
      DecompObjBound* objBoundLP = m_nodeStats.getLastBound();

      if (objBoundLP) {
         objBoundIP = *objBoundLP;
      }

      objBoundIP.thisBoundIP = thisBound;
      objBoundIP.bestBoundIP = m_nodeStats.objBest.second;
#ifdef UTIL_USE_TIMERS
      objBoundIP.timeStamp   = globalTimer.getRealTime();
#else
      objBoundIP.timeStamp   = -1;
#endif
      m_nodeStats.objHistoryBound.push_back(objBoundIP);
      UtilPrintFuncEnd(m_osLog, m_classTag,
                       "setObjBoundIP()", m_param.LogDebugLevel, 2);
   }


   bool isTailoffLB(const int    changeLen      = 10,
                    const double changePerLimit = 0.1);


   inline int getNumRowType(DecompRowType rowType) {
      int   nRowsType = 0;
      std::vector<DecompRowType>::iterator vi;

      for (vi = m_masterRowType.begin(); vi != m_masterRowType.end(); vi++) {
         if (*vi == rowType) {
            nRowsType++;
         }
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
              DecompApp*             app,
              UtilParameters*        utilParam):
      m_classTag   ("D-ALGO"),
      m_param      (),
      m_algo       (algo),
      m_status     (STAT_UNKNOWN),
      m_phase      (PHASE_UNKNOWN),
      m_phaseLast  (PHASE_UNKNOWN),
      m_phaseForce (PHASE_UNKNOWN),
      m_app        (app),
      m_stats      (),
      m_nodeStats  (),
      m_memPool    (),
      m_osLog      (&std::cout),
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
      m_masterObjLast(DecompInf),
      m_firstPhase2Call(false),
      m_isStrongBranch(false)
#ifdef DECOMP_MASTERONLY_DIRECT
      ,
      m_masterOnlyCols()
#endif
   {
      m_app->m_decompAlgo = this;
   }


   /**
    * Destructor.
    */
   virtual ~DecompAlgo() {
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
