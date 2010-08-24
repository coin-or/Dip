//[magala@orclus71 OsiCbc]$ g++ main.cpp -L ../../../build-g/lib -I ../../../build-g/include/coin/ -lCbc -lOsi -lCbc -lCgl -lCoinUtils -lOsiCbc -lClp -lOsi -lCoinUtils -lCbcSolver -lCbc -lOsiClp -lClp -lCgl -lCoinUtils -lOsi


#include "OsiCbcSolverInterface.hpp"
#include <string>
using namespace std;

int main(int argc, char ** argv){
   const int numCols      = 2;
   const int numNzs       = 10;
   const int numRows      =  6;
   bool      isRowOrdered = false;
   double    objective  [numCols] = {1,0};
   int       rowIndices [numNzs]  = {0,0,1,2,2,3,3,4,5,5};
   int       colIndices [numNzs]  = {0,1,1,0,1,0,1,1,0,1};
   double    elements   [numNzs]  = { 7.0, -1.0,  1.0, -1.0,  1.0,
                                      -4.0, -1.0, -1.0, 0.2, -1.0};
   CoinPackedMatrix M(isRowOrdered,
                      rowIndices, colIndices, elements, numNzs);
   double   rowLB[numRows] = {13.0, 1.0, -3.0, -27.0, -5.0, -4.0};
   double   rowUB[numRows] = {OsiCbcInfinity,
                              OsiCbcInfinity,
                              OsiCbcInfinity,
                              OsiCbcInfinity,
                              OsiCbcInfinity,
                              OsiCbcInfinity};
   double   colLB[numCols]       = {0,0};
   double   colUB[numCols]       = {6,6};
   int      integerVars[numCols] = {0,1};
   
   OsiCbcSolverInterface osi;
   osi.messageHandler()->setLogLevel(0);
   osi.loadProblem(M, colLB, colUB, objective, rowLB, rowUB);
   osi.setInteger(integerVars, 2);
   
   osi.branchAndBound();
   assert(!osi.isProvenPrimalInfeasible());
   assert(osi.isProvenOptimal());

   //osi-cbc changes internal column bounds, must reset
   double redCostX[numCols] = {1.34019,-0.10562};   
   osi.setColLower(colLB);
   osi.setColUpper(colUB);
   osi.setObjective(redCostX);
   osi.writeMps("tmp");
   osi.branchAndBound();
   assert(!osi.isProvenPrimalInfeasible());
   assert(osi.isProvenOptimal());

}
