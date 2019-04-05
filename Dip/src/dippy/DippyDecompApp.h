#ifndef DIPPY_DECOMPAPP_INCLUDED
#define DIPPY_DECOMPAPP_INCLUDED

//===========================================================================//
#include "Decomp.h"
#include "DecompApp.h"
#include "UtilParameters.h"

#include "Python.h"

#include <map>
#include <vector>
using namespace std;

//===========================================================================//
/*!
 * \class DippyDecompApp
 * A DecompApp that links Python to DIP.
 *
 * \see
 * DecompApp
 */

//===========================================================================//
class DippyDecompApp : public DecompApp {
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

   /** The various model constraint systems used for different algos. */

   PyObject* m_pProb;

   int m_numCols;

   /** Flags for Python callbacks. */
   bool m_pySolveRelaxed;
   bool m_pyIsUserFeasible;
   bool m_pyGenerateCuts;
   bool m_pyHeuristics;
   bool m_pyInitVars;

public:
   /** @name Helper functions (public). */

   /* Add PuLP problem. */

   void addPuLPProb(PyObject* p) {
      Py_XINCREF(p);          /* Add a reference to new problem */
      Py_XDECREF(m_pProb);    /* Dispose of previous problem */
      m_pProb = p;
   };

   /* Create models. */
   void createModels();

   virtual DecompSolverStatus solveRelaxed(const int whichBlock,
                                           const double* redCostX,
					   const double convexDual,
					   DecompVarList& varList);

   bool APPisUserFeasible(const double* x, const int n_cols, const double tolZero);

   virtual int generateCuts(const double* x, DecompCutList& newCuts);

   int APPheuristics(const double* xhat, const double* origCost, vector<DecompSolution*>& xhatIPFeas);

   int generateInitVars(DecompVarList& initVars);

   PyObject* m_rowList;
   map<PyObject*, int> m_rowIndices;
   PyObject* m_colList;
   map<PyObject*, int> m_colIndices;

   PyObject* m_relaxedKeys;
   map<PyObject*, int> m_relaxIndices;

public:
   DippyDecompApp(UtilParameters& utilParam, PyObject* p) :
      DecompApp  (utilParam),
      m_classTag ("SMALL-APP"),
      m_pProb(NULL) {
      addPuLPProb(p);
      createModels();
      m_pySolveRelaxed   = utilParam.GetSetting("pyRelaxedSolver", true);
      m_pyIsUserFeasible = utilParam.GetSetting("pyIsSolutionFeasible", true);
      m_pyGenerateCuts   = utilParam.GetSetting("pyGenerateCuts", true);
      m_pyHeuristics     = utilParam.GetSetting("pyHeuristics", true);
      m_pyInitVars       = utilParam.GetSetting("pyInitVars", true);
   }

   virtual ~DippyDecompApp() {
      // Remove references to Python objects/lists
      Py_XDECREF(m_pProb);
      Py_XDECREF(m_rowList);
      Py_XDECREF(m_colList);
      Py_XDECREF(m_relaxedKeys);
      delete [] m_objective;
      m_objective = NULL;
      delete m_modelCore.getModel();
      m_modelCore.setModel(NULL);
      map<int, DecompModel	>::iterator mit;

      for (mit = m_modelRelax.begin(); mit != m_modelRelax.end(); mit++) {
         DecompModel	& modelRelax = (*mit).second;
         delete modelRelax.getModel();
         modelRelax.setModel(NULL);
      }
   };
};

#endif
