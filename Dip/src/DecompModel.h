//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Authors: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)      //
//          Ted Ralphs, Lehigh University (ted@lehigh.edu)                   //
//          Jiadong Wang, Lehigh University (jiw408@lehigh.edu)              //
//                                                                           //
// Copyright (C) 2002-2015, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef DECOMP_MODEL_INCLUDED
#define DECOMP_MODEL_INCLUDED

//===========================================================================//
#include "UtilMacrosDecomp.h"
#include "DecompParam.h"
#include "DecompConstraintSet.h"
#include "DecompSolverResult.h"

//===========================================================================//
//naming convention - usually would do DecompModelXx, DecompModelYy
//===========================================================================//
class DecompModel {
protected:
   DecompConstraintSet*  m_model;
   std::string           m_modelName;
   int                   m_blockId;
   UtilParameters*       m_utilParam;

public:
   DecompConstraintSet* getModel()     const {
      return m_model;
   }
   const std::string&      getModelName() const {
      return m_modelName;
   }
   const int             getBlockId()   const {
      return m_blockId;
   }

public:
   void setModel    (DecompConstraintSet* model) {
      m_model = model;
   }
   void setModelName(const std::string modelName) {
      m_modelName = modelName;
   }
   void setBlockId  (const int blockId) {
      m_blockId = blockId;
   }

public:
   DecompModel(const DecompModel& appModel) {
      m_model     = appModel.m_model;
      m_modelName = appModel.m_modelName;
      m_blockId   = appModel.m_blockId;
      m_utilParam = appModel.m_utilParam;
   }
   
   DecompModel& operator=(const DecompModel& rhs) {
      m_model     = rhs.m_model;
      m_modelName = rhs.m_modelName;
      m_blockId   = rhs.m_blockId;
      m_utilParam = rhs.m_utilParam;
      return *this;
   }
   
   DecompModel(UtilParameters& utilParam) :
      m_model    (NULL),
      m_modelName(""),
      m_blockId  (0),
      m_utilParam(&utilParam){};

   DecompModel(DecompConstraintSet* model,
	       std::string          modelName,
	       int                  blockId,
	       UtilParameters&      utilParam) :
      m_model    (model),
      m_modelName(modelName),
      m_blockId  (blockId),
      m_utilParam(&utilParam){};

   virtual ~DecompModel() {}
};

//===========================================================================//
class DecompSubModel : public DecompModel {
private:
   OsiSolverInterface*   m_osi;
   int                   m_numCols;
   int*                  m_colIndices;
   int                   m_counter;
public:

   inline void setCounter(const int num) {
      m_counter = num;
   }

   inline int getCounter() {
      return m_counter;
   }

   void setOsi(OsiSolverInterface* osi) {
      m_osi = osi;

      if (!m_colIndices) {
         //---
         //--- For use with (re-)setting the objective coefficients
         //---  setup an array of indices 0,...,n-1. This object assumes
         //---  that the number of columns stayed constant throughout.
         //---
         const int numCols = m_osi->getNumCols();
         m_numCols         = numCols;
         m_colIndices      = new int[numCols];

         if (!m_colIndices) {
            UtilExceptionMemory("setOsi", "DecompSubModel");
         }

         UtilIotaN(m_colIndices, numCols, 0);
      }
   }

   void setOsiObjCoeff(const double* objCoeff) {
      assert(m_osi);
      assert(m_colIndices);
      assert(m_numCols == m_osi->getNumCols());

      if (getModel()->isSparse()) {
         const DecompConstraintSet* model = getModel();
         const std::map<int, int>& origToSparse = model->getMapOrigToSparse();
         std::map<int, int>::const_iterator mcit;

         for (mcit = origToSparse.begin();
               mcit != origToSparse.end(); mcit++) {
            m_osi->setObjCoeff(mcit->second,           //sparse-index
                               objCoeff[mcit->first]); //original-index
         }
      } else
         m_osi->setObjCoeffSet(m_colIndices,
                               m_colIndices + m_numCols, objCoeff);
   }


   void setActiveColBounds(const double* colLB,
                           const double* colUB) {
      DecompConstraintSet* model         = getModel();
      std::vector<int>&          activeColumns = model->activeColumns;

      //---
      //--- if no active columns are set,  assume they are all active
      //---   for e.g., in the case of one block (or sparse)
      //---
      if (model->isSparse()) {
         const std::map<int, int>& origToSparse = model->getMapOrigToSparse();
         std::map<int, int>::const_iterator mcit;

         for (mcit  = origToSparse.begin();
               mcit != origToSparse.end(); mcit++) {
            m_osi->setColLower(mcit->second,        //sparse-index
                               colLB[mcit->first]); //original-index
            m_osi->setColUpper(mcit->second,        //sparse-index
                               colUB[mcit->first]); //original-index
            //printf("setColBounds orig:%d sparse:%d lb:%g ub:%g\n",
            //   mcit->first, mcit->second,
            //   colLB[mcit->first],
            //   colUB[mcit->first]);
         }
      } else {
         if (activeColumns.size()) {
            std::vector<int>::iterator vi;

            for (vi = activeColumns.begin(); vi != activeColumns.end(); vi++) {
               //printf("setColBounds i:%d to LB:%g UB:%g\n",
               //     *vi, colLB[*vi], colUB[*vi]);
               m_osi->setColBounds(*vi, colLB[*vi], colUB[*vi]);
            }
         } else {
            m_osi->setColLower(colLB);
            m_osi->setColUpper(colUB);
         }
      }
   }

   void solveAsMIPSym(DecompSolverResult*  result,
		      DecompParam&         param,
		      bool                 doExact,
		      bool                 doCutoff,
		      bool                 isRoot,
		      double               cutoff,
		      double               timeLimit);
   
   void solveAsMIPCbc(DecompSolverResult*  result,
		      DecompParam&         param,
		      bool                 doExact,
		      bool                 doCutoff,
		      bool                 isRoot,
		      double               cutoff,
		      double               timeLimit);
   
   void solveAsMIPCpx(DecompSolverResult*  result,
		      DecompParam&         param,
		      bool                 doExact,
		      bool                 doCutoff,
		      bool                 isRoot,
		      double               cutoff,
		      double               timeLimit);
   
   void solveAsMIPGrb(DecompSolverResult*  result,
		      DecompParam&         param,
		      bool                 doExact,
		      bool                 doCutoff,
		      bool                 isRoot,
		      double               cutoff,
		      double               timeLimit);
   
public:
   OsiSolverInterface*   getOsi() const {
      return m_osi;
   }

public:
   void solveAsMIP(DecompSolverResult*  result,
		   DecompParam&         param,
		   bool                 doExact,
		   bool                 doCutoff,
		   bool                 isRoot,
		   double               cutoff,
		   double               timeLimit);
   
   bool isPointFeasible(const double*  x,
                        const bool     isXSparse  = false,
                        const int      logLevel   = 0,
                        const double   feasVarTol = 1.0e-5,
                        const double   feasConTol = 1.0e-4);

public:
   DecompSubModel(const DecompModel& appModel) :
      DecompModel(appModel),
      m_osi         (NULL),
      m_numCols     (0   ),
      m_colIndices  (NULL),
      m_counter     ( 0 )
   {};

   DecompSubModel& operator=(const DecompModel& rhs) {
      DecompModel::operator=(rhs);
      return *this;
   }

   DecompSubModel(UtilParameters &utilParam) :
      DecompModel(utilParam),
      m_osi         (NULL),
      m_numCols     (0   ),
      m_colIndices  (NULL),
      m_counter     (0)
   {};
   DecompSubModel(DecompConstraintSet* model,
                  std::string          modelName,
		  int                  blockId,
		  UtilParameters&      utilParam) :
      DecompModel(model, modelName, blockId, utilParam),
      m_osi         (NULL),
      m_numCols     (0   ),
      m_colIndices  (NULL),
      m_counter     (0)
   {};
   ~DecompSubModel() {
      if (m_osi) {
         delete    m_osi;
      }

      if (m_colIndices) {
         delete [] m_colIndices;
      }
   }
};

#endif
