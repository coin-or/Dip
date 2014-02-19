#ifndef DIPPY_DECOMPALGO_INCLUDED
#define DIPPY_DECOMPALGO_INCLUDED

#include "Decomp.h"
#include "DecompAlgoC.h"
#include "DecompAlgoPC.h"
#include "DecompAlgoRC.h"
#include "DecompAlgoD.h"
#include "DecompCutPool.h"
#include "AlpsDecompTreeNode.h"

#include "Python.h"

/**
 * Mixin class for Dip Algorithms
 *
 * This is a helper class for interfacing Dip Algo classes with Python. To add Python support
 * to a standard DecompAlgo, create a subclass which also inherits from DippyAlgoMixin and
 * override the virtual methods to call those provided by the Mixin class. See DippyAlgoC
 * for an example.
 */
class DippyAlgoMixin {
private:
   PyObject* pDownLB;
   PyObject* pDownUB;
   PyObject* pUpLB;
   PyObject* pUpUB;

protected:
   PyObject* m_pProb;
   UtilParameters* m_utilParam;
public:
   /**
    * Constructor
    *
    * @param utilParam parameter class
    * @param pProb a DipProblem python object
    */
   DippyAlgoMixin(UtilParameters* utilParam, PyObject* pProb)
      : m_pProb(pProb), m_utilParam(utilParam), pDownLB(NULL),
        pDownUB(NULL), pUpLB(NULL), pUpUB(NULL) {
   }

   bool chooseBranchSet(DecompAlgo* algo,
                        std::vector< std::pair<int, double> >& downBranchLB,
                        std::vector< std::pair<int, double> >& downBranchUB,
                        std::vector< std::pair<int, double> >& upBranchLB,
                        std::vector< std::pair<int, double> >& upBranchUB);

   PyObject* getPDownLB() {
      return pDownLB;
   };
   PyObject* getPDownUB() {
      return pDownUB;
   };
   PyObject* getPUpLB() {
      return pUpLB;
   };
   PyObject* getPpUpUB() {
      return pUpUB;
   };

   void postProcessBranch(DecompAlgo* algo, DecompStatus decompStatus);

   void postProcessNode(DecompAlgo* algo, DecompStatus decompStatus);

   /**
    * Generate Cuts
    *
    * Call the algorithm's generateCuts method if the 'generateCuts' param is set to
    * true.
    */
   virtual int generateCuts(DecompAlgo* algo, double* xhat, DecompCutList& newCuts) {
      bool gen = m_utilParam->GetSetting("generateCuts", true);
      bool inPy = m_utilParam->GetSetting("pyGenerateCuts", true);

      if (gen && inPy) {
         return algo->DecompAlgo::generateCuts(xhat, newCuts);
      }

      return 0;
   }

   /**
    * Generate Initial Variables
    *
    * Call the algorithm's generateInitVars method if the 'generateInitVars' param is set to
    * true.
    */
   virtual int generateInitVars(DecompAlgo* algo, DecompVarList& initVars) {
      bool gen = m_utilParam->GetSetting("generateInitVars", true);
      bool inPy = m_utilParam->GetSetting("pyInitVars", true);

      if (gen && inPy) {
         return algo->DecompAlgo::generateInitVars(initVars);
      }

      return 0;
   }
};

/**
 * Python-enabled DecompAlgoC
 *
 */
class DippyAlgoC : DippyAlgoMixin, public DecompAlgoC {
public:
   DippyAlgoC(DecompApp* app, UtilParameters* utilParam, PyObject* pProb)
      : DippyAlgoMixin(utilParam, pProb), DecompAlgoC(app, utilParam) {
   }

   virtual bool chooseBranchSet(
      std::vector< std::pair<int, double> >& downBranchLB,
      std::vector< std::pair<int, double> >& downBranchUB,
      std::vector< std::pair<int, double> >& upBranchLB,
      std::vector< std::pair<int, double> >& upBranchUB) {
      bool ret_val = DippyAlgoMixin::chooseBranchSet(this, downBranchLB, downBranchUB,
                     upBranchLB, upBranchUB);
      return ret_val;
   }

   virtual void postProcessBranch(DecompStatus decompStatus) {
      DippyAlgoMixin::postProcessBranch(this, decompStatus);
   }

   virtual void postProcessNode(DecompStatus decompStatus) {
      DippyAlgoMixin::postProcessNode(this, decompStatus);
   }

   virtual int generateCuts(double* xhat, DecompCutList& newCuts) {
      return DippyAlgoMixin::generateCuts(this, xhat, newCuts);
   }

   //    virtual int generateInitVars(DecompVarList & initVars) {
   //	      return priceInitVars(initVars);
   //    }
};

/**
 * Python-enabled DecompAlgoPC
 *
 */
class DippyAlgoPC : DippyAlgoMixin, public DecompAlgoPC {
public:
   DippyAlgoPC(DecompApp* app, UtilParameters* utilParam, PyObject* pProb)
      : DippyAlgoMixin(utilParam, pProb), DecompAlgoPC(app, utilParam) {
   }

   virtual bool chooseBranchSet(std::vector< std::pair<int, double> >& downBranchLB,
                                std::vector< std::pair<int, double> >& downBranchUB,
                                std::vector< std::pair<int, double> >& upBranchLB,
                                std::vector< std::pair<int, double> >& upBranchUB) {
      return DippyAlgoMixin::chooseBranchSet(this, downBranchLB, downBranchUB,
                                             upBranchLB, upBranchUB);
   }

   virtual void postProcessBranch(DecompStatus decompStatus) {
      DippyAlgoMixin::postProcessBranch(this, decompStatus);
   }

   virtual void postProcessNode(DecompStatus decompStatus) {
      DippyAlgoMixin::postProcessNode(this, decompStatus);
   }

   virtual int generateCuts(double* xhat, DecompCutList& newCuts) {
      return DippyAlgoMixin::generateCuts(this, xhat, newCuts);
   }

   virtual int generateInitVars(DecompVarList& initVars) {
      return DippyAlgoMixin::generateInitVars(this, initVars);
   }
};

/**
 * Python-enabled DecompAlgoRC
 *
 */
class DippyAlgoRC : DippyAlgoMixin, public DecompAlgoRC {
public:
   DippyAlgoRC(DecompApp* app, UtilParameters* utilParam, PyObject* pProb)
      : DippyAlgoMixin(utilParam, pProb), DecompAlgoRC(app, utilParam) {
   }

   virtual bool chooseBranchSet(std::vector< std::pair<int, double> >& downBranchLB,
                                std::vector< std::pair<int, double> >& downBranchUB,
                                std::vector< std::pair<int, double> >& upBranchLB,
                                std::vector< std::pair<int, double> >& upBranchUB) {
      return DippyAlgoMixin::chooseBranchSet(this, downBranchLB, downBranchUB,
                                             upBranchLB, upBranchUB);
   }

   virtual void postProcessBranch(DecompStatus decompStatus) {
      DippyAlgoMixin::postProcessBranch(this, decompStatus);
   }

   virtual void postProcessNode(DecompStatus decompStatus) {
      DippyAlgoMixin::postProcessNode(this, decompStatus);
   }

   virtual int generateCuts(double* xhat, DecompCutList& newCuts) {
      return DippyAlgoMixin::generateCuts(this, xhat, newCuts);
   }

   virtual int generateInitVars(DecompVarList& initVars) {
      return DippyAlgoMixin::generateInitVars(this, initVars);
   }
};

#endif
