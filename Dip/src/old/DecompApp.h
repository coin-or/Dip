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


#ifndef DECOMP_APP_INCLUDED
#define DECOMP_APP_INCLUDED

/*-----------------------------------------------------------------------*/
#include "UtilParameters.h"
#include "DecompParam.h"
#include "DecompModel.h"
#include "DecompSolution.h"
#include "DecompConstraintSet.h"

/*-----------------------------------------------------------------------*/
/*!
 * \class DecompApp
 * \brief
 * The main application class.
 *
 * The main application class where the user will define the model
 * decomposition and define any application specific methods.
 *
 * \see
 * DecompModel
 * DecompConstraintSet
 * DecompParam
 *
 * \todo
 * Clone a monkey.
 */
/*-----------------------------------------------------------------------*/
class DecompApp {

private:
   /*!
    * Disable copy constructors.
    */
   DecompApp(const DecompApp&);
   DecompApp& operator=(const DecompApp&);

private:
   /*!
    * Store the name of the class (for logging/debugging) - "who am I?"
    */
   static const char* m_classTag;

protected:
   /*!
    *  Log file.
    */
   ostream* m_osLog;

public:
   /*!
    *  Parameters.
    */
   DecompParam m_param;

   /*!
    *  Model data object.
    */
   DecompModel                    m_model;
   map<int, DecompConstraintSet*> m_modelCore;
   map<int, DecompConstraintSet*> m_modelRelax;

public:
   //set log file
   /**
      Initialize the DecompApp data.
   */
   void startupLog();

   //base layer needs to do some kind of check to make sure this actually
   //got done - but also nice to have user version... so base createModel
   //and user uesrCreateModel() which is pure?, in base userCreateModel
   //gets called and checked that it returns good information
   int createModel();
   //interface method idea of Gardner - so it is very clear
   //what they can and cannot override and what they must
   virtual void APPcreateModel(double                        *& objCoeff,
                               map<int, DecompConstraintSet*> & modelCore,
                               map<int, DecompConstraintSet*> & modelRelax) = 0;

   //?
   virtual bool APPisUserFeasible(const double* x,
                                  const int      n_cols,
                                  const double   tolZero) {
      return true;
   };

   virtual int APPheuristics(const double*              xhat,
                             vector<DecompSolution*>  & xhatIPFeas) {
      return 1;
   }

   virtual int generateInitVars(DecompVarList& initVars,
                                int whichModel);

   //APPgenerateCuts?
   virtual int generateCuts(const double*               x,
                            const DecompConstraintSet& modelCore,
                            const DecompConstraintSet& modelRelax,
                            DecompCutList&              newCuts);

   virtual decompStat APPsolveRelaxed(const int             whichModel,
                                      const double*         redCostX,
                                      const double*         origCost,
                                      const double          alpha,
                                      const int             n_origCols,
                                      const bool            checkRC,
                                      const bool            checkDup,
                                      OsiSolverInterface*   m_subprobSI,
                                      list<DecompVar*>    & vars) {
      return STAT_UNKNOWN;
   }


   virtual void printOriginalColumn(const int   index,
                                    ostream*    os = &cout) const;

   virtual void printOriginalSolution(const int      n_cols,
                                      const double* solution,
                                      ostream* os = &cout) const;

   //Cut(this){
public:
   DecompApp(UtilParameters& utilParam) :
      m_osLog(&cout),
      m_param(),
      m_model(),
      m_modelCore(),
      m_modelRelax() {
      m_param.getSettings(utilParam);
      startupLog();
   };
   virtual ~DecompApp() {
      map<int, DecompConstraintSet*>::iterator it;

      for (it = m_modelCore.begin(); it != m_modelCore.end(); it++) {
         UTIL_DELPTR(it->second);
      }

      for (it = m_modelRelax.begin(); it != m_modelRelax.end(); it++) {
         UTIL_DELPTR(it->second);
      }
   };
};

#endif
