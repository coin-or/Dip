
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


#ifndef DECOMP_APP_INCLUDED
#define DECOMP_APP_INCLUDED

//===========================================================================//
#include "UtilParameters.h"
#include "DecompParam.h"
#include "DecompModel.h"
#include "DecompSolution.h"
#include "DecompConstraintSet.h"
#include "CoinMpsIO.hpp"
#include "CoinLpIO.hpp"

extern "C" {
#if defined (COIN_HAS_METIS)
#include "hmetis.h"
#endif
}
//===========================================================================//
class DecompAlgo;

//===========================================================================//
/*!
 * \class DecompApp
 * \brief
 * The main application class.
 *
 * The main application class where the user will define the model
 * decomposition and define any application specific methods.
 *
 */

//===========================================================================//
class DecompApp {

private:
   /**
    * Store the name of the class (for logging/debugging) - "who am I?"
    */
   std::string m_classTag;

protected:
   /**
    *  Log file.
    */
   std::ostream* m_osLog;

   /**
    * The best known LB/UB for this application (if known, for debugging).
    */
   double m_bestKnownLB;
   double m_bestKnownUB;





public:
   /**
    *  Parameters.
    */
   DecompParam m_param;

   /**
    *  Model data: objective function
    */
   const double* m_objective;


   /**
    *  Model data: the core model (A'')
    */
   DecompAppModel m_modelCore;

   /**
    *  Model data: the relaxed model(s) (A')
    */
   std::map<int, DecompAppModel> m_modelRelax;

   /**
    *  Model data: the relaxed (nested) model(s) (A')
    */
   std::map<int, std::vector<DecompAppModel> > m_modelRelaxNest;

   /**
    * Pointer to the base algorithmic object.
    *   NOTE: only for the advanced user
    */
   DecompAlgo* m_decompAlgo;


   /*
    *  The following definitiions was from MILPBlock
    *
    */

   /** MPS object for reading instances */
   CoinMpsIO m_mpsIO;

   /** LP object for reading instances */
   CoinLpIO m_lpIO;

   /** Original constraint matrix for the instance */

   const CoinPackedMatrix* m_matrix;


   /** The model constraint systems used for different algos */

   DecompConstraintSet*   m_modelC;
   std::map<int, DecompConstraintSet*>  m_modelR;


   /** Definition of blocks (by rows) */

   std::map<int, std::vector<int> > m_blocks;



public:
   /**
    * @name Helper functions.
    * @{
    */


   /** Preprocess (standard ): on the TODO list*/

   void preprocess();

   /**
    * Print startup message to log.
    */
   void startupLog();

   //TODO:
   //base layer needs to do some kind of check to make sure this actually
   //got done - but also nice to have user version... so base createModel
   //and user uesrCreateModel() which is pure?, in base userCreateModel
   //gets called and checked that it returns good information
   int createModel();

   inline const double getBestKnownLB() const {
      return m_bestKnownLB;
   }
   inline const double getBestKnownUB() const {
      return m_bestKnownUB;
   }

   inline void setBestKnownLB(const double bestKnownLB) {
      m_bestKnownLB = bestKnownLB;
   }
   inline void setBestKnownUB(const double bestKnownUB) {
      m_bestKnownUB = bestKnownUB;
   }


   /**
    * Set the model objective function.
    *
    * NOTE: The user application MUST call this method.
    */
   inline void setModelObjective(const double* objective) {
      assert(objective);
      m_objective = objective;
   }

   /**
    * Set the model core constraint matrix.
    *
    * NOTE: The user application MUST call this method.
    */
   //TODO: having these implementations in the header makes
   // it harder to view this as an interface class - it is unclear
   // what the user must do vs can do
   inline void setModelCore(DecompConstraintSet* model,
                            const std::string          modelName) {
      assert(model);

      if (!model->hasPrepRun()) {
         model->prepareModel(true);
      }

      m_modelCore.setModel(model);
      m_modelCore.setModelName(modelName);
   }
   ////////////STOP thoughts on redesign
   //TODO: change to setModelCore( ... )
   //  to long set of args of basic types - more like how Osi
   //  does loadProblem, etc... or juse use Osi as shell and accept
   //  an OsiSolverInterface?? or an OsiModel... is there one?
   //   maybe use CoinModel?? much cleaner?
   //at least set up an alternative to do it that way
   //still, this is still a function the user must call versus
   // a method they MUST override -which would fit more into the
   // framework concept - for e.g.,
   //virtual CoinModel * createModelCore() - the use must create a CoinModel
   //  or maybe a DecompModel which is derived from a CoinModel with
   //  whatever extra stuff might be needed - but I would try to avoid that
   //is it best to have it returned as return of function - forcing the
   //user to do it? or as arguments? and copy or assign pointers like
   //in SAS load problem design



   /**
    * Set the model relaxed constraint matrix (for a particular block).
    *
    * NOTE: The user application MUST call this method IF they are
    * not deriving the function DecompApp::solveRelaxed.
    */
   inline void setModelRelax(DecompConstraintSet* model,
                             const std::string          modelName = "",
                             const int             blockId   = 0) {
      if (model && !model->hasPrepRun()) {
         model->prepareModel();
      }

      //---
      //--- make sure this block has not been set yet
      //---
      std::map<int, DecompAppModel>::iterator mit = m_modelRelax.find(blockId);

      if (mit != m_modelRelax.end()) {
         std::cerr << "Block " << blockId << " relaxation has already been set. "
                   << "Only one relaxation definition can be used at one time."
                   << std::endl;
         throw UtilException("Multiple relaxation definitions",
                             "setModelRelax", "DecompApp");
      }

      DecompAppModel appModel(model, modelName, blockId);
      m_modelRelax.insert(std::make_pair(blockId, appModel));
   }

   /**
    * Set the model relaxed (nested) constraint matrix
    *  (for a particular block).
    */
   inline void setModelRelaxNest(DecompConstraintSet* model,
                                 const std::string          modelName,
                                 const int             blockId = 0) {
      assert(model);

      if (!model->hasPrepRun()) {
         model->prepareModel();
      }

      DecompAppModel appModel(model, modelName, blockId);
      m_modelRelaxNest[blockId].push_back(appModel);
   }

   /**
    * Get a pointer to the base algorithm class
    */
   inline DecompAlgo* getDecompAlgo() const  {
      return m_decompAlgo;
   }


   /**
    * @}
    */


public:
   /**
    * @name Interface methods for user derivation (virtual).
    * @{
    */

   /**
    * Initialize the dual vector for PhaseII of PC. The user is passed
    * a reference to the internal data and can manipulate it directly.
    *
    * This is only called when dual stabilization is used, i.e.,
    * when m_param.DualStab > 0, at the first iteration of PhaseII of PC.
    * The vector is immediately smoothed with the initial restricted master
    * duals. By default, the restricted mater is used as the initial dual
    * and, therefore, no smoothing occurs in the first iteration.
   */
   virtual void initDualVector(std::vector<double>& dualVector) {}


   //TODO: change name - no other one is using APP, why here?
   /**
    * Method to determine if the solution (x) is feasible to the original
    * model. For explicitly defined model components, like the model core
    * constraints (A''), the feasibility of the solution is automatically
    * checked against the constraints. In the case when the relaxed problem
    * constraints (A') are explicitly defined - these are also checked
    * automatically.
    *
    * However, for some applications, a valid feasible constraint system
    * cannot be explicitly defined (even for the core set of constraints).
    * For example, think of the case of TSP, where A'' is defined as the
    * subtour elimination constraints. These constraints are implicitly
    * defined by deriving the method DecompApp::generateCuts. Therefore,
    * the framework cannot automatically tell if a solution is feasible
    * by checking against the constraint system. In this case, the user must
    * provide this method.
    *
    * @param[in] x The solution point to check.
    * @return True, if x is feasible; otherwise, false.
    */
   //TODO: what is doxy tag for function return
   //TOOD: don't need numCols and tolZero should not be user overidable
   virtual bool APPisUserFeasible(const double* x,
                                  const int      numCols,
                                  const double   tolZero) {
      return true;
   };

   virtual int APPheuristics(const double*             xhat,
                             const double*             origCost,
                             std::vector<DecompSolution*>& xhatIPFeas) {
      return 0;
   }

   /**
    * This function allows the user to return their own dual vector
    * to be used in the generation of new variables (in the reduced-cost
    * calculation).
    *
    * For reference, the user is given the dual vector from the restricted
    * master (or the stabilized dual, if using m_param.DualStab).
    */
   virtual const double* getDualForGenerateVars(const double* dual) {
      return 0;
   }

   virtual int generateInitVars(DecompVarList& initVars);

   virtual int generateCuts(const double*   x,
                            DecompCutList& newCuts);

   virtual void solveRelaxedWhich(std::vector<int>&                blocksToSolve,
                                  std::map< int, std::vector<double> >& userDualsByBlock) {
   };

   virtual DecompSolverStatus solveRelaxed(const int          whichBlock,
                                           const double*      redCostX,
                                           DecompVarList&     varList) {
      return DecompSolStatNoSolution;
   }
   virtual DecompSolverStatus solveRelaxedNest(const int          whichBlock,
         const double*      redCostX,
         DecompVarList&     varList) {
      return DecompSolStatNoSolution;
   }


   virtual void printOriginalColumn(const int   index,
                                    std::ostream*    os = &std::cout) const;

   //TODO: change api so colNames comes from modelCore if exists
   //  rather than - to simplify API
   virtual void printOriginalSolution(const int              n_cols,
                                      const std::vector<std::string>& colNames,
                                      const double*          solution,
                                      std::ostream*               os = &std::cout) const;

   /**
    * @}
    */



public:

   /** Initialize applications */
   virtual void initializeApp(UtilParameters& utilParam);

   /** Create model parts */

   void createModels();

   DecompConstraintSet* createModelPart(const int nRowsPart,
                                        const int* rowsPart);

   void createModelPart(DecompConstraintSet* model,
                        const int             nRowsPart,
                        const int*            rowsPart);

   void createModelPartSparse(DecompConstraintSet* model,
                              const int             nRowsPart,
                              const int*            rowsPart);

   void createModelMasterOnlys(std::vector<int>& masterOnlyCols);

   void readInitSolutionFile(DecompVarList& initVars);

   /** Read block file */
   void readBlockFile();

   /** Read Problem */
   const CoinPackedMatrix*  readProblem(UtilParameters& utilParam);


   /** Automatically detect singly bordered structure */

   void singlyBorderStructureDetection();

   /** Find the active columns for some block */
   void findActiveColumns(const std::vector<int>& rowsPart,
                          std::set<int>&           activeColsSet);

   /** Get Intance name */

   const std::string getInstanceName() {
      return m_param.Instance;
   }
   /*
   void  HMETIS_PartKway(int nvtxs, int nhedges, int *vwgts, int *eptr,
   	      int *eind, int *hewgts, int nparts, int ubfactor,
   	      int * options, int * part, int *edgecut);

   void  HMETIS_PartRecursive(int nvtxs, int nhedges, int *vwgts, int *eptr,
   		   int *eind, int *hewgts, int nparts, int ubfactor,
   		   int * options, int * part, int *edgecut);
   */




public:
   /**
    * Constructor for base DecompApp class. This accepts a generic parameters
    * object (UtilParameters) and reads in the parameter settings into the
    * DecompApp paramter object.
    */
   DecompApp() :
      m_classTag   ("D-APP"),
      m_osLog      (&std::cout  ),
      m_bestKnownLB(-1e75  ),
      m_bestKnownUB( 1e75  ),
      m_objective  ( 0     ),
      m_matrix     ( 0     ),
      m_modelC     ( 0     ),
      m_modelR     ( ) {
      //---
      //--- comment these functions, which were used in
      //--- MILPBlock, otherwise, conflict occurs in building
      //--- individual examples
      //     m_param.getSettings(utilParam);
      //     initializeApp(utilParam);
      //     startupLog();
   };

   /**
    * Destructor.
    */
   virtual ~DecompApp() {
      UTIL_DELARR(m_objective);
      UtilDeleteMapPtr(m_modelR);
      UTIL_DELPTR(m_modelC);
      UTIL_DELPTR(m_matrix);
      //     std::cout << "destructor is called " << std::endl;
   };
};

#endif
