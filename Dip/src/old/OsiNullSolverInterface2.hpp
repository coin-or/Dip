// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiNullSolverInterface_H
#define OsiNullSolverInterface_H

#include <string>
#include <vector>

#include "CoinMessageHandler.hpp"
#include "CoinPackedVectorBase.hpp"

#include "OsiData.hpp"
#include "OsiCollections.hpp"
#include "OsiSolverParameters.hpp"

class CoinPackedMatrix;
class CoinWarmStart;

class OsiCuts;
class OsiAuxInfo;
class OsiRowCut;
class OsiRowCutDebugger;
class CoinSet;
class CoinBuild;
class CoinModel;
class OsiSolverBranch;
class OsiSolverResult;
#include "CoinFinite.hpp"

#ifndef COIN_DBL_MAX
static const double OsiNullInfinity = DBL_MAX;
#else
static const double OsiNullInfinity = COIN_DBL_MAX;
#endif

//#############################################################################

class OsiNullSolverInterface : public OsiSolverInterface  {

public:
   ///@name Solve methods 
   //@{
   /// Solve initial LP relaxation 
   virtual void initialSolve(){
      CoinAssertHint(0, "OsiNull does not have a solver");
   }; 
   
   /// Resolve an LP relaxation after problem modification
   virtual void resolve(){
      CoinAssertHint(0, "OsiNull does not have a solver");
   }
   
   /// Invoke solver's built-in enumeration algorithm
   virtual void branchAndBound(){
      CoinAssertHint(0, "OsiNull does not have a solver");
   }
   //@{
   // Set an integer parameter
   // copy all parameters in this section from one solver to another
   void copyParameters(OsiNullSolverInterface & rhs);
   //@}

   //------------------------------------------------------------------------
   ///@name Methods returning info on how the solution process terminated
   //@{
   /// Are there numerical difficulties?
   virtual bool isAbandoned() const{
      CoinAssertHint(0, "OsiNull does not have a solver");
      return false;
   }
   /// Is optimality proven?
   virtual bool isProvenOptimal() const{
      CoinAssertHint(0, "OsiNull does not have a solver");
      return false;
   }
   /// Is primal infeasiblity proven?
   virtual bool isProvenPrimalInfeasible() const{
      CoinAssertHint(0, "OsiNull does not have a solver");
      return false;
   }
   /// Is dual infeasiblity proven?
   virtual bool isProvenDualInfeasible() const{
      CoinAssertHint(0, "OsiNull does not have a solver");
      return false;
   }
   /// Is the given primal objective limit reached?
   virtual bool isPrimalObjectiveLimitReached() const{
      CoinAssertHint(0, "OsiNull does not have a solver");
      return false;
   }
   /// Is the given dual objective limit reached?
   virtual bool isDualObjectiveLimitReached() const{
      CoinAssertHint(0, "OsiNull does not have a solver");
      return false;
   }
   /// Iteration limit reached?
   virtual bool isIterationLimitReached() const{
      CoinAssertHint(0, "OsiNull does not have a solver");
      return false;
   }
   //@}

   //------------------------------------------------------------------------
   /** \name Warm start methods

   Note that the warm start methods return a generic CoinWarmStart object.
   The precise characteristics of this object are solver-dependent. Clients
   who wish to maintain a maximum degree of solver independence should take
   care to avoid unnecessary assumptions about the properties of a warm start
   object.
   */
   //@{
   /*! \brief Get an empty warm start object
      
   This routine returns an empty warm start object. Its purpose is
   to provide a way for a client to acquire a warm start object of the
   appropriate type for the solver, which can then be resized and modified
   as desired.
   */

   virtual CoinWarmStart *getEmptyWarmStart () const {
      CoinAssertHint(0, "OsiNull does not have a solver");
      return NULL;
   }

   /** \brief Get warm start information.

   Return warm start information for the current state of the solver
   interface. If there is no valid warm start information, an empty warm
   start object wil be returned.
   */
   virtual CoinWarmStart* getWarmStart() const {
      CoinAssertHint(0, "OsiNull does not have a solver");
      return NULL;
   }

   /** \brief Set warm start information.
    
   Return true or false depending on whether the warm start information was
   accepted or not.
   By definition, a call to setWarmStart with an empty warm start object
   should remove the warm start information held in the solver interface.
   */
   virtual bool setWarmStart(const CoinWarmStart* warmstart){
      CoinAssertHint(0, "OsiNull does not have a solver");
      return false;
   }
   //@}

   //------------------------------------------------------------------------
   /**@name Problem query methods

   Querying a problem that has no data associated with it will result in
   zeros for the number of rows and columns, and NULL pointers from
   the methods that return vectors.
     
   Const pointers returned from any data-query method are valid as
   long as the data is unchanged and the solver is not called.
   */
   //@{
   /// Get number of columns
   virtual int getNumCols() const {
      return data_->getNcol();
   };
  
   /// Get number of rows
   virtual int getNumRows() const {
      return data_->getNrow();
   };
  
   /// Get number of nonzero elements
   virtual int getNumElements() const {
      const CoinPackedMatrix * Mrow = data_->getMatrixByRow();
      const CoinPackedMatrix * Mcol = data_->getMatrixByCol();
      if(Mrow){
         return Mrow->getNumElements();
      } else if(Mcol){
         return Mcol->getNumElements();
      }
      return 0;
   };

   /// Get number of integer variables
   virtual int getNumIntegers() const {
      int i, nInt = 0;
      const char * colType = data_->getColType();
      for(i = 0; i < getNumCols(); i++){
         if(colType[i] == 'B' || colType[i] == 'I')
            nInt++;
      }
      return nInt;
   }
  
   /// Get pointer to array[getNumCols()] of column lower bounds
   virtual const double * getColLower() const {
      return data_->getColLower();
   }
   
   /// Get pointer to array[getNumCols()] of column upper bounds
   virtual const double * getColUpper() const {
      return data_->getColUpper();
   }
   
   /** Get pointer to array[getNumRows()] of row constraint senses.
       <ul>
       <li>'L': <= constraint
       <li>'E': =  constraint
       <li>'G': >= constraint
       <li>'R': ranged constraint
       <li>'N': free constraint
       </ul>
   */
   virtual const char * getRowSense() const {
  
      /** Get pointer to array[getNumRows()] of row right-hand sides
          <ul>
          <li> if getRowSense()[i] == 'L' then
          getRightHandSide()[i] == getRowUpper()[i]
          <li> if getRowSense()[i] == 'G' then
          getRightHandSide()[i] == getRowLower()[i]
          <li> if getRowSense()[i] == 'R' then
          getRightHandSide()[i] == getRowUpper()[i]
          <li> if getRowSense()[i] == 'N' then
          getRightHandSide()[i] == 0.0
          </ul>
      */
      return data_->getRowSense();
   }

   virtual const double * getRightHandSide() const {
  
      /** Get pointer to array[getNumRows()] of row ranges.
          <ul>
          <li> if getRowSense()[i] == 'R' then
          getRowRange()[i] == getRowUpper()[i] - getRowLower()[i]
          <li> if getRowSense()[i] != 'R' then
          getRowRange()[i] is 0.0
          </ul>
      */
      return data_->getRowRhs();
   }

   virtual const double * getRowRange() const {
      return data_->getRowRange();
   }

   /// Get pointer to array[getNumRows()] of row lower bounds
   virtual const double * getRowLower() const {
      return data_->getRowLower();
   }
  
   /// Get pointer to array[getNumRows()] of row upper bounds
   virtual const double * getRowUpper() const {
      return data_->getRowUpper();
   }
  
   /// Get pointer to array[getNumCols()] of objective function coefficients
   virtual const double * getObjCoefficients() const {
      return data_->getObj();
   }
  
   /// Get objective function sense (1 for min (default), -1 for max)
   virtual double getObjSense() const {
      CoinAssertHint(0, "Sorry, not implemented yet.");
      return 1;
   }

   /// Return true if variable is continuous
   virtual bool isContinuous(int colIndex) const {
      return data_->getColType()[colIndex] == 'C';
   }
  
   /// Return true if variable is binary
   virtual bool isBinary(int colIndex) const {
      return data_->getColType()[colIndex] == 'B';
   }
  
   /** Return true if column is integer.
       Note: This function returns true if the the column
       is binary or a general integer.
   */
   virtual bool isInteger(int colIndex) const {
      return (data_->getColType()[colIndex] == 'B' || 
              data_->getColType()[colIndex] == 'I' ); 
   }
    
   /// Get pointer to row-wise copy of matrix
   virtual const CoinPackedMatrix * getMatrixByRow() const {
      return data_->getMatrixByRow();
   }
  
   /// Get pointer to column-wise copy of matrix
   virtual const CoinPackedMatrix * getMatrixByCol() const {
      return data_->getMatrixByCol();
   }
   /// Get solver's value for infinity
   virtual double getInfinity() const {
      return data_->getInfinity();
   }
   //@}
    
   /**@name Solution query methods */
   //@{
   /// Get pointer to array[getNumCols()] of primal variable values
   virtual const double * getColSolution() const {
      return data_->getPrimalSol();
   }
  
   /// Get pointer to array[getNumRows()] of dual variable values
   //TODO: let user enter dual sol also
   virtual const double * getRowPrice() const {
      CoinAssertHint(0, "Sorry, not implemented yet.");
      return NULL;
   }
  
   /// Get a pointer to array[getNumCols()] of reduced costs
   virtual const double * getReducedCost() const {
      CoinAssertHint(0, "Sorry, not implemented yet.");
      return NULL;
   }
  
   /** Get pointer to array[getNumRows()] of row activity levels (constraint
       matrix times the solution vector). */
   virtual const double * getRowActivity() const {
      return data_->getRowActivity();
   }
  
   /// Get objective function value
   virtual double getObjValue() const {
      CoinAssertHint(0, "Sorry, not implemented yet.");
      return 0.0;
   }
   
   /** Get the number of iterations it took to solve the problem (whatever
       ``iteration'' means to the solver). */
   virtual int getIterationCount() const {
      CoinAssertHint(0, "OsiNull does not have a solver");
   }
  
   /** Get as many dual rays as the solver can provide. In case of proven
       primal infeasibility there should be at least one.
     
       \note
       Implementors of solver interfaces note that
       the double pointers in the vector should point to arrays of length
       getNumRows() and they should be allocated via new[].
     
       \note
       Clients of solver interfaces note that
       it is the client's responsibility to free the double pointers in the
       vector using delete[].
   */
   virtual std::vector<double*> getDualRays(int maxNumRays) const {
      /** Get as many primal rays as the solver can provide. (In case of proven
          dual infeasibility there should be at least one.)
     
          <strong>NOTE for implementers of solver interfaces:</strong> <br>
          The double pointers in the vector should point to arrays of length
          getNumCols() and they should be allocated via new[]. <br>
     
          <strong>NOTE for users of solver interfaces:</strong> <br>
          It is the user's responsibility to free the double pointers in the
          vector using delete[].
      */
      CoinAssertHint(0, "OsiNull does not have a solver");
   }

   virtual std::vector<double*> getPrimalRays(int maxNumRays) const {
      CoinAssertHint(0, "OsiNull does not have a solver");
   }
  
   //-------------------------------------------------------------------------
   /**@name Methods to modify the objective, bounds, and solution

   For functions which take a set of indices as parameters
   (\c setObjCoeffSet(), \c setColSetBounds(), \c setRowSetBounds(),
   \c setRowSetTypes()), the parameters follow the C++ STL iterator
   convention: \c indexFirst points to the first index in the
   set, and \c indexLast points to a position one past the last index
   in the set.    
   */
   //@{

   /** Set an objective function coefficient */
   virtual void setObjCoeff( int elementIndex, double elementValue ) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Set a single column lower bound.
       Use -getInfinity() for -infinity. */
   virtual void setColLower( int elementIndex, double elementValue ) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Set a single column upper bound.
       Use getInfinity() for infinity. */
   virtual void setColUpper( int elementIndex, double elementValue ) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }
  
   /** Set a single row lower bound.
       Use -getInfinity() for -infinity. */
   virtual void setRowLower( int elementIndex, double elementValue ) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }
      
   /** Set a single row upper bound.
       Use getInfinity() for infinity. */
   virtual void setRowUpper( int elementIndex, double elementValue ) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }
    
   /** Set the type of a single row */
   virtual void setRowType(int index, char sense, double rightHandSide,
                           double range) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /// Set the objective function sense.
   /// (1 for min (default), -1 for max)
   virtual void setObjSense(double s) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Set characters for columns types.
       colType[i] can have values
       <UL>
       <LI> 'C' : continuous
       <LI> 'B' : binary
       <LI> 'I' : integer
       </UL>
   */
   virtual void setColType(const char * colType){
      data_->setColType(colType);
   }
  
   /*
     Set the primal solution variable values
    
     colsol[getNumCols()] is an array of values for the primal variables.
     These values are copied to memory owned by the solver interface object
     or the solver.  They will be returned as the result of getColSolution()
     until changed by another call to setColSolution() or by a call to any
     solver routine.  Whether the solver makes use of the solution in any
     way is solver-dependent.
   */
   virtual void setColSolution(const double *colsol) {
      data_->setPrimalSol(colsol);
   }


   /** Set dual solution variable values

   rowprice[getNumRows()] is an array of values for the dual
   variables. These values are copied to memory owned by the solver
   interface object or the solver.  They will be returned as the result of
   getRowPrice() until changed by another call to setRowPrice() or by a
   call to any solver routine.  Whether the solver makes use of the
   solution in any way is solver-dependent.
   */

   virtual void setRowPrice(const double * rowprice) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   //-------------------------------------------------------------------------
   /**@name Methods to set variable type */
   //@{
   /** Set the index-th variable to be a continuous variable */
   virtual void setContinuous(int index) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Set the index-th variable to be an integer variable */
   virtual void setInteger(int index) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Set the variables listed in indices (which is of length len) to be
       continuous variables */
   //@}
   //-------------------------------------------------------------------------
    
   //-------------------------------------------------------------------------
   /**@name Methods to modify the constraint system.

   Note that new columns are added as continuous variables.

   */
   //@{
   /** Add a column (primal variable) to the problem. */
   virtual void addCol(const CoinPackedVectorBase& vec,
                       const double collb, const double colub,   
                       const double obj) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** \brief Remove a set of columns (primal variables) from the
       problem.

       The solver interface for a basis-oriented solver will maintain valid
       warm start information if all deleted variables are nonbasic.
   */
   virtual void deleteCols(const int num, const int * colIndices) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }
    
   /** Add a row (constraint) to the problem. */
   virtual void addRow(const CoinPackedVectorBase& vec,
                       const double rowlb, const double rowub) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** */
   virtual void addRow(const CoinPackedVectorBase& vec,
                       const char rowsen, const double rowrhs,   
                       const double rowrng) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** \brief Delete a set of rows (constraints) from the problem.

   The solver interface for a basis-oriented solver will maintain valid
   warm start information if all deleted rows are loose.
   */
   virtual void deleteRows(const int num, const int * rowIndices) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }
    


   //---------------------------------------------------------------------
   /**@name Methods to input a problem */
   //TODO: ugh force both row/col format?
   void loadDataAndSolution(const CoinPackedMatrix & rowMatrix,
                            const CoinPackedMatrix & colMatrix,
                            const double           * collb, 
                            const double           * colub,   
                            const double           * obj,
                            const double           * rowlb, 
                            const double           * rowub,
                            const char             * colType,
                            const double           * primalSol,
                            const double             infinity){
      data_->setMatrixByRow(&rowMatrix);
      data_->setMatrixByCol(&colMatrix);      
      data_->setNrow(rowMatrix.getNumRows());
      data_->setNcol(rowMatrix.getNumCols());
      data_->setColLower(collb);
      data_->setColUpper(colub);
      data_->setRowLower(rowlb);
      data_->setRowUpper(rowub);
      data_->setObj(obj);
      data_->setColType(colType);
      data_->setPrimalSol(primalSol);
      data_->setInfinity(infinity);
      data_->initializeOtherData();
   }



   //@{
   /** Load in an problem by copying the arguments (the constraints on the
       rows are given by lower and upper bounds). If a pointer is 0 then the
       following values are the default:
       <ul>
       <li> <code>colub</code>: all columns have upper bound infinity
       <li> <code>collb</code>: all columns have lower bound 0 
       <li> <code>rowub</code>: all rows have upper bound infinity
       <li> <code>rowlb</code>: all rows have lower bound -infinity
       <li> <code>obj</code>: all variables have 0 objective coefficient
       </ul>
   */
   virtual void loadProblem(const CoinPackedMatrix & matrix,
                            const double * collb, 
                            const double * colub,   
                            const double * obj,
                            const double * rowlb, 
                            const double * rowub){
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   //TODO: version that passes in rhs/ranges/sense and converts to bounds
			    
   /** Load in an problem by assuming ownership of the arguments (the
       constraints on the rows are given by lower and upper bounds).
       For default values see the previous method.

       \warning
       The arguments passed to this method will be
       freed using the C++ <code>delete</code> and <code>delete[]</code>
       functions. 
   */
   virtual void assignProblem(CoinPackedMatrix*& matrix,
                              double*& collb, double*& colub, double*& obj,
                              double*& rowlb, double*& rowub) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Load in an problem by copying the arguments (the constraints on the
       rows are given by sense/rhs/range triplets). If a pointer is 0 then the
       following values are the default:
       <ul>
       <li> <code>colub</code>: all columns have upper bound infinity
       <li> <code>collb</code>: all columns have lower bound 0 
       <li> <code>obj</code>: all variables have 0 objective coefficient
       <li> <code>rowsen</code>: all rows are >=
       <li> <code>rowrhs</code>: all right hand sides are 0
       <li> <code>rowrng</code>: 0 for the ranged rows
       </ul>
   */
   virtual void loadProblem(const CoinPackedMatrix& matrix,
                            const double* collb, const double* colub,
                            const double* obj,
                            const char* rowsen, const double* rowrhs,   
                            const double* rowrng) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Load in an problem by assuming ownership of the arguments (the
       constraints on the rows are given by sense/rhs/range triplets). For
       default values see the previous method.

       \warning
       The arguments passed to this method will be
       freed using the C++ <code>delete</code> and <code>delete[]</code>
       functions. 
   */
   virtual void assignProblem(CoinPackedMatrix*& matrix,
                              double*& collb, double*& colub, double*& obj,
                              char*& rowsen, double*& rowrhs,
                              double*& rowrng) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Just like the other loadProblem() methods except that the matrix is
       given in a standard column major ordered format (without gaps). */
   virtual void loadProblem(const int numcols, const int numrows,
                            const CoinBigIndex * start, const int* index,
                            const double* value,
                            const double* collb, const double* colub,   
                            const double* obj,
                            const double* rowlb, const double* rowub) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Just like the other loadProblem() methods except that the matrix is
       given in a standard column major ordered format (without gaps). */
   virtual void loadProblem(const int numcols, const int numrows,
                            const CoinBigIndex * start, const int* index,
                            const double* value,
                            const double* collb, const double* colub,   
                            const double* obj,
                            const char* rowsen, const double* rowrhs,   
                            const double* rowrng) {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }

   /** Write the problem in MPS format to the specified file.

   If objSense is non-zero, a value of -1.0 causes the problem to be
   written with a maximization objective; +1.0 forces a minimization
   objective. If objSense is zero, the choice is left to implementation.
   */
   virtual void writeMps(const char *filename,
                         const char *extension = "mps",
                         double objSense=0.0) const {
      CoinAssertHint(0, "Sorry, not implemented yet.");
   }
   
   //-----------------------------------------------------------------------

   ///@name Constructors and destructors
   //@{
   /// Default Constructor
   OsiNullSolverInterface() :
      OsiSolverInterface(),
      data_(NULL)
   {
      data_ = new OsiData(OsiNullInfinity);
   }
   
   /** Clone

   The result of calling clone(false) is defined to be equivalent to
   calling the default constructor OsiNullSolverInterface().
   */
   virtual OsiNullSolverInterface * clone(bool copyData = true) const {
      CoinAssertHint(0, "Sorry, not implemented yet.");
      return NULL;
   }
  
   /// Copy constructor (disabled)
   OsiNullSolverInterface(const OsiNullSolverInterface &);
  
   /// Assignment operator (disabled)
   OsiNullSolverInterface & operator=(const OsiNullSolverInterface& rhs);
  
   /// Destructor 
   virtual ~OsiNullSolverInterface () {
      delete data_;
   }

   //@}

   //------------------------------------------------------------------------

protected:
   ///@name Protected methods
   //@{
   /** Apply a row cut (append to the constraint matrix). */
   virtual void applyRowCut( const OsiRowCut & rc ) {
      CoinAssertHint(0, "OsiNull does not have a solver");
   }

   /** Apply a column cut (adjust the bounds of one or more variables). */
   virtual void applyColCut( const OsiColCut & cc ) {
      CoinAssertHint(0, "OsiNull does not have a solver");
   }
   /** A quick inlined function to force a value to be between a minimum and
       a maximum value */
   template <class T> inline T
   forceIntoRange(const T value, const T lower, const T upper) const {
      return value < lower ? lower : (value > upper ? upper : value);
   }
   //@}
  
   //---------------------------------------------------------------------
protected:
   OsiData * data_;
};

#endif
