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
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
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
class DecompApp{

private:
   /**
    * Store the name of the class (for logging/debugging) - "who am I?"
    */
   string m_classTag;
  
protected:
   /**
    *  Log file.
    */
   ostream * m_osLog;

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
   const double * m_objective;

   /**
    *  Model data: the core model (A'')
    */  
   DecompAppModel m_modelCore;

   /**
    *  Model data: the relaxed model(s) (A')
    */  
   map<int, DecompAppModel> m_modelRelax;

   /**
    *  Model data: the relaxed (nested) model(s) (A')
    */  
   map<int, vector<DecompAppModel> > m_modelRelaxNest;   

   /**
    * Pointer to the base algorithmic object.
    *   NOTE: only for the advanced user
    */
   DecompAlgo * m_decompAlgo;
   
public:
   /**
    * @name Helper functions.
    * @{
    */

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

   inline const double getBestKnownLB() const {return m_bestKnownLB;}
   inline const double getBestKnownUB() const {return m_bestKnownUB;}

   inline void setBestKnownLB(const double bestKnownLB) {
      m_bestKnownLB = bestKnownLB; 
   }
   inline void setBestKnownUB(const double bestKnownUB) {
      m_bestKnownUB = bestKnownUB; 
   }


   //should run some sanity checks here!
   inline void setModelObjective(const double * objective){
      assert(objective);
      m_objective = objective;
   }
   inline void setModelCore(DecompConstraintSet * model,
                            const string          modelName){
      assert(model);
      if(!model->hasPrepRun())
         model->prepareModel();
      m_modelCore.setModel(model);
      m_modelCore.setModelName(modelName);
   }

   inline void setModelRelax(DecompConstraintSet * model,
                             const string          modelName = "",
                             const int             blockId   = 0){
      if(model && !model->hasPrepRun())
	 model->prepareModel();      
      
      //---
      //--- make sure this block has not been set yet
      //---
      map<int, DecompAppModel>::iterator mit = m_modelRelax.find(blockId);
      if(mit != m_modelRelax.end()){
	 cerr << "Block " << blockId << " relaxation has already been set. "
	      << "Only one relaxation definition can be used at one time." 
	      << endl;
	 throw UtilException("Multiple relaxation definitions",
			     "setModelRelax", "DecompApp");
      }	 
      DecompAppModel appModel(model, modelName, blockId);
      m_modelRelax.insert(make_pair(blockId, appModel));
   }
   
   inline void setModelRelaxNest(DecompConstraintSet * model,
                                 const string          modelName,
                                 const int             blockId = 0){
      //assuming blocks are disjoint in variables - if not, bug
      //   can check that with active columns
      assert(model);
      if(!model->hasPrepRun())
         model->prepareModel();
      DecompAppModel appModel(model, modelName, blockId);
      m_modelRelaxNest[blockId].push_back(appModel);
   }
   
   /**
    * Get a pointer to the base algorithm class
    */
   inline DecompAlgo * getDecompAlgo() const  {
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
      Initialize the dual vector for PhaseII of PC. 
      
      This is only called when dual stabilization is used, i.e.,
      when m_param.DualStab > 0.      
   */
   //virtual void initDualVector(vector<double> & dualVector){
   //}


   virtual bool APPisUserFeasible(const double * x,
                                  const int      n_cols,
                                  const double   tolZero){
      return true;
   };

   virtual int APPheuristics(const double            * xhat,
			     const double            * origCost,
			     vector<DecompSolution*> & xhatIPFeas){
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
   virtual const double * getDualForGenerateVars(const double * dual){
      return 0;
   }
  
   virtual int generateInitVars(DecompVarList & initVars);
   
   virtual int generateCuts(const double  * x, 
			    DecompCutList & newCuts);

   virtual void solveRelaxedWhich(vector<int> & blocksToSolve){};
      
   virtual DecompSolverStatus solveRelaxed(const int          whichBlock,
                                           const double     * redCostX,
                                           const double       convexDual,
                                           DecompVarList    & varList){
      return DecompSolStatNoSolution;
   }
   virtual DecompSolverStatus solveRelaxedNest(const int          whichBlock,
                                               const double     * redCostX,
                                               const double       convexDual,
                                               DecompVarList    & varList){
      return DecompSolStatNoSolution;
   }


   virtual void printOriginalColumn(const int   index, 
                                    ostream   * os = &cout) const;
  
   //TODO: change api so colNames comes from modelCore if exists
   //  rather than - to simplify API 
   virtual void printOriginalSolution(const int              n_cols, 
                                      const vector<string> & colNames,
                                      const double         * solution,
                                      ostream              * os = &cout) const;

   /**
    * @}
    */

   

public:
   /**
    * Constructor for base DecompApp class. This accepts a generic parameters
    * object (UtilParameters) and reads in the parameter settings into the 
    * DecompApp paramter object. 
    */
   DecompApp(UtilParameters & utilParam) :
      m_classTag   ("D-APP"),
      m_osLog      (&cout  ), 
      m_bestKnownLB(-1e75  ),
      m_bestKnownUB( 1e75  ),
      m_objective  (0      ),
      m_decompAlgo (0      )
   {    
      m_param.getSettings(utilParam);
      startupLog();    
   };

   /**
    * Destructor.
    */
   virtual ~DecompApp() {};
};

#endif
