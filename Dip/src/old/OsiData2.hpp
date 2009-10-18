// Name:     OsiData.hpp
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
//           email: fmargot@andrew.cmu.edu
// Date:     12/2/06
//-----------------------------------------------------------------------------
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.

#ifndef OsiData_H
#define OsiData_H

/** Class collecting pointers on data for OsiEmpty. Each generator
    may have a derived class to add additional pointers on data. 
    If a data member is not used by
    a generator, the data member need not be defined (or may be NULL).
    Ownership of the data remains with the calling method.
*/

#include "CoinPackedMatrix.hpp"

class OsiData {

public:

   /**@name Public Set/get methods */
   //@{

   /** Set infinity. */
   virtual void setInfinity(const double givenInfinity){
      infinity = givenInfinity;
   }
   
   /** Get infinity */
   inline double getInfinity() const {return infinity;};

   /** Set nrow to the number of rows */
   virtual void setNrow(const int givenNrow){
      nrow = givenNrow;
   }
   
   /** Get nrow */
   inline int getNrow() const {return nrow;};

   /** Set ncol to the number of variables */
   virtual void setNcol(const int givenNcol){
      ncol = givenNcol;
   }

   /** Get ncol */
   inline int getNcol() const {return ncol;};

   /** Set  matrixByCol to point on the coefficient matrix ordered by 
       columns */
   virtual void setMatrixByCol(const CoinPackedMatrix *givenMatrixByCol){
        matrixByCol = givenMatrixByCol;
   }

   /** Get matrixByCol */
   inline const CoinPackedMatrix *getMatrixByCol() const {return matrixByCol;};

   /** Set  matrixByRow to point on the coefficient matrix ordered by 
       rows */
   virtual void setMatrixByRow(const CoinPackedMatrix *givenMatrixByRow){
      matrixByRow = givenMatrixByRow;
   }

   /** Get matrixByRow */
   inline const CoinPackedMatrix *getMatrixByRow() const {return matrixByRow;};

   /** Set obj to point on a vector holding the objective coefficient values */
   virtual void setObj(const double *givenObj){
      obj = givenObj;
   }

   /** Get obj */
   inline const double *getObj() const {return obj;};

   /** Set  colLower to point on a vector holding the lower bounds on the 
       variables */
   virtual void setColLower(const double *givenColLower){
      colLower = givenColLower;
   }

   /** Get colLower */
   inline const double *getColLower() const {return colLower;};

   /** Set  colUpper to point on a vector holding the upper bounds on the 
       variables */
   virtual void setColUpper(const double *givenColUpper){
      colUpper = givenColUpper;
   }
   
   /** Get colUpper */
   inline const double *getColUpper() const {return colUpper;};

   /** Set  rowLower to point on a vector holding the lower bounds on the 
       constraints */
   virtual void setRowLower(const double *givenRowLower){
      rowLower = givenRowLower;
   }
   
   /** Get rowLower */
   inline const double *getRowLower() const {return rowLower;};

   /** Set  rowUpper to point on a vector holding the upper bounds on the 
       constraints */
   virtual void setRowUpper(const double *givenRowUpper){
      rowUpper = givenRowUpper;
   }

   /** Get rowUpper */
   inline const double *getRowUpper() const {return rowUpper;};

   /** Set  rowRhs to point on a vector holding the right hand side of the 
       constraints (for a ranged constraint, it contains the upper bound). */
   //virtual void setRowRhs(const double *givenRowRhs);
   /** Get rowRhs */
   inline const double *getRowRhs() const {return rowRhs;};
   /** Get rowRange */
   inline const double *getRowRange() const {return rowRange;};
   /** Get rowSense */
   inline const char *getRowSense() const {return rowSense;};

   /** Set  rowActivity to point on a vector holding the activity of the 
       constraints (i.e. coefficient matrix times separateThis). */
   //virtual void setRowActivity(const double *givenRowActivity);
   /** Get rowActivity */
   inline const double *getRowActivity() const {return rowActivity;};

   /** Set colType to point on a vector holding the type of the
       variables ('B', 'I', or 'C' for Binary, Integer and Continuous) */ 
   virtual void setColType(const char *givenColType){
      colType = givenColType;
   }

   /** Get colType */
   inline const char *getColType() const {return colType;};

   /** Set primal solution. */
   virtual void setPrimalSol(const double * givenPrimalSol){
      primalSol = givenPrimalSol;
   }

   /** Get primal solution. */
   inline const double *getPrimalSol() const {return primalSol;};

   /** initialize the non-const data */
   void initializeOtherData(){
      CoinAssert(matrixByRow);
      CoinAssert(primalSol);
      
      rowRhs      = new double[nrow];
      rowRange    = new double[nrow];
      rowSense    = new char[nrow];
      rowActivity = new double[nrow];
      
      int    r;
      for(r = 0; r < nrow; r++){
         convertBoundToSense(rowLower[r], rowUpper[r], 
                             rowSense[r], rowRhs[r], rowRange[r]);
      }
      matrixByRow->times(primalSol, rowActivity);
   }
   //@}

   /** A quick inlined function to convert from the lb/ub style of constraint
       definition to the sense/rhs/range style */
   inline void convertBoundToSense(const double   lower, 
                                   const double   upper,
                                   char         & sense, 
                                   double       & right,
                                   double       & range) const
   {
      double inf = getInfinity();
      range = 0.0;
      if (lower > -inf) {
         if (upper < inf) {
            right = upper;
            if (upper==lower) {
               sense = 'E';
            } else {
               sense = 'R';
               range = upper - lower;
            }
         } else {
            sense = 'G';
            right = lower;
         }
      } else {
         if (upper < inf) {
            sense = 'L';
            right = upper;
         } else {
            sense = 'N';
            right = 0.0;
         }
      }
   }

   /**@name Constructors and destructors */
   //@{
   /// Default constructor 
   OsiData(const double givenInfinity = DBL_MAX,
           const int &givenNrow = 0, 
           const int &givenNcol = 0,
           const CoinPackedMatrix * givenMatrixByCol = NULL,
           const CoinPackedMatrix * givenMatrixByRow = NULL,
           const double * givenObj = NULL,
           const double * givenColLower = NULL, 
           const double * givenColUpper = NULL,   
           const double * givenRowLower = NULL, 
           const double * givenRowUpper = NULL,
           const char   * givenColType = NULL,
           const double * givenPrimalSol = NULL) :
      infinity(givenInfinity),
      nrow(givenNrow),
      ncol(givenNcol),
      matrixByCol(givenMatrixByCol),
      matrixByRow(givenMatrixByRow),
      obj(givenObj),
      colLower(givenColLower),
      colUpper(givenColUpper),
      rowLower(givenRowLower),
      rowUpper(givenRowUpper),
      colType(givenColType),
      primalSol(givenPrimalSol) {}
   
   
   /// Destructor 
   virtual ~OsiData(){
      if(rowRhs)
         delete [] rowRhs;
      if(rowRange)
         delete [] rowRange;
      if(rowSense)
         delete [] rowSense;
      if(rowActivity)
         delete [] rowActivity;
   }
      //@}

protected:

   // Private member data

   /**@name Private member data */

   //@{

   // Definition of infinity
   double infinity;

   // Number of constraints
   int nrow;

   // Number of variables.
   int ncol;

   // Pointer on matrix of coefficients (ordered by columns).
   CoinPackedMatrix const *matrixByCol;

   // Pointer on matrix of coefficients (ordered by rows).
   CoinPackedMatrix const *matrixByRow;

   // Pointer on vector of objective coefficients. 
   const double *obj;

   // Pointer on vector of lower bounds on variables.
   const double *colLower; 

   // Pointer on vector of upper bounds for variables.
   const double *colUpper;   

   // Pointer on vector of lower bounds for constraints.
   const double *rowLower; 

   // Pointer on vector of upper bounds for constraints.
   const double *rowUpper;

   /** Pointer on vector of characters for columns types.
       colType[i] can have values
       <UL>
       <LI> 'C' : continuous
       <LI> 'B' : binary
       <LI> 'I' : integer
       </UL>
   */
   const char * colType;

   // Pointer on vector for current primal solution.
   const double * primalSol;

protected:
   
   // Vector of right hand sides for constraints.
   double * rowRhs;

   // Vector of ranges for constraints.
   double * rowRange;

   // Vector of row senses for constraints.
   char   * rowSense;

   // Vector of activity of constraints (i.e. coefficient matrix 
   // times primalSol).
   double * rowActivity;

   
   //@}

};



#endif
