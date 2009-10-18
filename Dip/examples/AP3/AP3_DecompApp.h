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

#ifndef AP3_DECOMPAPP_INCLUDED
#define AP3_DECOMPAPP_INCLUDED

// --------------------------------------------------------------------- //
#include "DecompApp.h"
#include "AP3_Instance.h"
#include "AP3_DecompParam.h"
// --------------------------------------------------------------------- //

/*!
 * \class AP3_DecompApp
 * A DecompApp for solving the 3-Indexed Assignment Problem (AP3).
 * 
 * \see
 * DecompApp
 *
 */

// --------------------------------------------------------------------- //
class AP3_DecompApp : public DecompApp{
private:
   /** Class id tag (for log / debugging). */
   const string m_classTag;

private:
   /** AP3 problem instance data */
   AP3_Instance m_ap3data;
   
private:
   /** Helper data for fast cut separation. */
   map<int, vector<int> > m_intersection;

   /** Helper data for fast construction of relaxed problems. */
   OsiSolverInterface  * m_siAP;
   double             ** m_assigncostMin;
   int                ** m_assignindexMin;
   
protected:
   /** Application specific parameters. */
   AP3_DecompParam m_appParam;  
   
public:
   /** Model types */
   enum ModelType {
      /** Relaxation is 2-D Assignment Problem over JxK */
      MODEL_I,
      /** Relaxation is 2-D Assignment Problem over IxK */
      MODEL_J,
      /** Relaxation is 2-D Assignment Problem over IxJ */
      MODEL_K
   };
   
public:
   /* @name Inherited (from pure virtual) methods. */

   /** Create the application model(s). */
   //TODO: model object? vs objCoeff?
   void APPcreateModel(double                        *& objCoeff,
                       map<int, DecompConstraintSet*> & modelCore,
                       map<int, vector<DecompConstraintSet*> > & modelRelax);

public:
   /* @name Inherited (from virtual) methods. */

   /** Solve the relaxed problem. */
   //TOOD: too messy?
   DecompStatus APPsolveRelaxed(const int             whichModel,
                                const double        * redCostX,
                                const double        * origCost,
                                const double          alpha,
                                const int             n_origCols,
                                const bool            checkRC,
                                const bool            checkDup,
                                OsiSolverInterface  * m_subprobSI,
                                list<DecompVar*>    & vars);
   
   /** Print an original column. */
   void printOriginalColumn(const int   index,
                            ostream   * os = &cout) const ;

public:
   /** @name Helper functions (public). */

   /** Get optimal (or best known) objective. */
   inline const double getKnownOptimalBound() const {
     return m_ap3data.m_optBound;
   }

   /** Get name of instance. */
   inline const string getInstanceName() const {
     return m_ap3data.m_instance;
   }
   
   /** Global index for column (i,j,k) in 3D case. */
   inline int indexIJK(const int i, 
                       const int j, 
                       const int k) const {
      int n = m_ap3data.m_dimension;
      return (i * n * n) + (j * n) + k;
   }
   inline int indexJIK(const int j, 
                       const int i, 
                       const int k) const {
      int n = m_ap3data.m_dimension;
      return (i * n * n) + (j * n) + k;
   }
   inline int indexKIJ(const int k, 
                       const int i, 
                       const int j) const {
      int n = m_ap3data.m_dimension;
      return (i * n * n) + (j * n) + k;
   }

   /** Global index for column (ind1, ind2, ind3) in 3D case. */
   inline void index3Inv(const int   index,
                         int       & ind1, 
                         int       & ind2,
                         int       & ind3) const { 
      int n   = m_ap3data.m_dimension;
      int nsq = n * n;
      int indexMod;
      ind1     = index / nsq;
      indexMod = index % nsq; 
      ind2     = indexMod / n;
      ind3     = indexMod % n;
   }

   /** Global index for column (ind1, ind2) in 2D case. */
   inline int index2(const int ind1, 
                     const int ind2) const { 
      int n = m_ap3data.m_dimension;
      return (ind1 * n) + ind2;
   }

   /** Return the indices (ind1, ind2) for 2D case. */
   inline pair<int,int> index2Inv(const int index) const { 
      int n = m_ap3data.m_dimension;
      return make_pair(index / n, index % n);
   }

   /** Return the known optimal bound. */
   inline const double getOptBound() const{
      return m_ap3data.m_optBound;
   }
   
private: 
   /** @name Helper functions (private). */   

   /** Guts of constructor. */
   void initializeApp(UtilParameters & utilParam) throw(CoinError);

   /** Create constraint system for one relaxation. */
   void createModelPart(const int             modelType,
                        int                 * rowInd,
                        double              * rowEls,
                        DecompConstraintSet * modelCore,
                        DecompConstraintSet * modelRelax) throw(CoinError);
   
private:
   /** @name Copy Constructors
    *
    * Disable the default copy constructors.
    *
    */
   AP3_DecompApp(const AP3_DecompApp &);
   AP3_DecompApp & operator=(const AP3_DecompApp &);

public:
   /** @name Constructor and Destructor */

   /** Default constructor. Takes an instance of UtilParameters */
   AP3_DecompApp(UtilParameters & utilParam) : 
      DecompApp(utilParam),
      m_classTag("AP3-APP"),
      m_ap3data(),
      m_siAP(0),
      m_assigncostMin(0),
      m_assignindexMin(0),
      m_appParam()
   {
      initializeApp(utilParam);         
   }
   
   virtual ~AP3_DecompApp() {
      int d;
      for(d = 0; d < m_ap3data.m_dimension; d++){
         UTIL_DELARR(m_assigncostMin[d]);
         UTIL_DELARR(m_assignindexMin[d]);
      }
      UTIL_DELPTR(m_siAP);
      UTIL_DELARR(m_assigncostMin);
      UTIL_DELARR(m_assignindexMin);
   };
};

#endif
