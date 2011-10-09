// [magala@orclus71 bugs]$ g++ main.cpp -L ../../build-g/lib -I ../../build-g/incli -lClp -lCbc -lCoinUtils -lCgl -lOsi -lCbcSolver -lCbc -lCgl -lClp


#include "CbcSolver.hpp"
#include "OsiClpSolverInterface.hpp"
#include <string>
using namespace std;

int main(int argc, char ** argv){
   string lpFile = argv[1];

   OsiClpSolverInterface si;
   si.readLp(lpFile.c_str());

   CbcModel cbc(si);
   CbcMain0(cbc);

   const  char * cbcArgv[20];
   int    cbcArgc   = 0;
   string cbcExe    = "cbc";
   string cbcSolve  = "-solve";
   string cbcQuit   = "-quit";
   string cbcLog    = "-log";
   string cbcLogSet = "3";
   cbcArgv[cbcArgc++]  = cbcExe.c_str();
   cbcArgv[cbcArgc++]  = cbcLog.c_str();
   cbcArgv[cbcArgc++]  = cbcLogSet.c_str();      
   cbcArgv[cbcArgc++]  = cbcSolve.c_str();
   cbcArgv[cbcArgc++]  = cbcQuit.c_str();
   CbcMain1(cbcArgc, cbcArgv, cbc);
   printf("cbc.status() = %d\n", cbc.status());
   printf("cbc.isProveOptimal() = %d\n", cbc.isProvenOptimal());
}
