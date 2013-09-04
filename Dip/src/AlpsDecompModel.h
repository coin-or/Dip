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

//===========================================================================//
#ifndef AlpsDecompModel_h_
#define AlpsDecompModel_h_

//===========================================================================//
#include "Alps.h"
#include "AlpsModel.h"
#include "AlpsDecompParam.h"


//===========================================================================//
#include "DecompAlgo.h"
#include "DecompConstraintSet.h"

//===========================================================================//
#include <string>

//===========================================================================//
class AlpsTreeNode;

//===========================================================================//
/**
 * \class AlpsDecompModel
 * \brief Derivation of AlpsModel for DECOMP.
 *
 * An object derived from AlpsModel. It interfaces with DECOMP methods
 * through a pointer to the active DecompAlgo.
 *
 * - AlpsDecompModel is derived from AlpsModel
 *    - AlpsModel has no pure virtual functions
 * - AlpsModel is derived from AlpsKnowledge
 *    - AlpsKnowledge has no pure virtual functions
 *
 * Virtual methods that should be derived here:
 * - createRoot
 *
 * \see AlpsModel
 * \see DecompAlgo
 *
 * \todo Clone a monkey.
 * \todo Arm wrestle Ted.
 * \todo Allow use of Alps writeParameters.
 * \todo Use message handler.
 * \todo Use differencing scheme.
 * \todo Setup for parallel.
 */
//===========================================================================//

//===========================================================================//
class AlpsDecompModel : public AlpsModel {


private:

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
    * Parameters for Alps.
    */
   AlpsDecompParam m_param;

   /**
    * Pointer to decomp algorithm used for bounding.
    */
   DecompAlgo* m_decompAlgo;

   /**
    * Objective of best solution found.
    */
   double m_bestLB;
   double m_bestUB;
   int    m_nodesProcessed;
   int    m_alpsStatus;

   /**
    * @}
    */

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
   AlpsDecompModel() :
      AlpsModel    (),
      m_classTag   ("ALPSM"),
      m_param      (),
      m_decompAlgo (NULL) {
   }

   AlpsDecompModel(UtilParameters& utilParam,
                   DecompAlgo*      decompAlgo) :
      AlpsModel   (),
      m_classTag  ("ALPSM"),
      m_param     (utilParam),
      m_decompAlgo(decompAlgo) {
      if (decompAlgo == NULL)
         throw UtilException("No DecompAlgo algorithm has been set.",
                             "AlpsDecompModel", "AlpsDecompModel");

      setAlpsSettings();
   }

   /**
    * Destructor.
    */
   virtual ~AlpsDecompModel() {}
   /**
    * @}
    */


   //-----------------------------------------------------------------------//
   /**
    * @name Virtual functions from AlpsModel.
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * Create the root node of the search tree.
    */
   virtual AlpsTreeNode* createRoot();

   /** Return true, if all nodes can be fathomed.*/
   virtual bool fathomAllNodes();
   /**
    * @}
    */

   //-----------------------------------------------------------------------//
   /**
    * @name Helper functions.
    * @{
    */
   //-----------------------------------------------------------------------//
   /**
    * Solve with ALPS and DECOMP.
    */
   AlpsExitStatus solve();

   /**
    * Set the ALPS parameters.
    */
   void setAlpsSettings();

   /**
    * Solve with ALPS and DECOMP.
    */
   inline void setDecompAlgo(DecompAlgo* decompAlgo) {
      m_decompAlgo = decompAlgo;
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
   /**
    * Get a ptr to the decomp algorithm vector.
    */
   //THINK: return ref?
   inline DecompAlgo* getDecompAlgo() {
      return m_decompAlgo;
   }

   inline AlpsDecompParam& getParam() {
      return m_param;
   }


   /**
    * Get number of rows in core decomp model.
    */
   inline const int getNumCoreRows() const {
      return m_decompAlgo->getModelCore().getModel()->getNumRows();
   }

   /**
    * Get number of cols in core decomp model.
    */
   inline const int getNumCoreCols() const {
      return m_decompAlgo->getModelCore().getModel()->getNumCols();
   }

   /**
    * Get the column names in core decomp model.
    */
   inline const std::vector<std::string>& getColNames() const {
      return m_decompAlgo->getModelCore().getModel()->getColNames();
   }

   /**
    * Get the row names in core decomp model.
    */
   inline const std::vector<std::string>& getRowNames() const {
      return m_decompAlgo->getModelCore().getModel()->getRowNames();
   }

   /**
    * Get the best solution found.
    */
   inline const DecompSolution* getBestSolution() const {
      return m_decompAlgo->getXhatIPBest();
   }

   const double getGlobalLB() const {
      return m_bestLB;
   }
   const double getGlobalUB() const {
      return m_bestUB;
   }
   const int    getSolStatus() const {
      return m_alpsStatus;
   }
   const int    getNumNodesProcessed() const {
      return m_nodesProcessed;
   }

   /**
    * @}
    */
};

#endif
