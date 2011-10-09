//#if 0
//--------------------------------------------------------------------- //
//for debugging
bool MMKP_DecompApp::APPisUserFeasible(const double * x,
				       const int      n_cols,
				       const double   tolZero){

   //---
   //--- Assume: it is already integeral.
   //---  s.t. sum{i in 1..n, j in 1..l[i]}  r[k,i,j] x[i,j] <= b[k], k in 1..m
   //---       sum{j in 1..l[i]}                      x[i,j]  = 1   , i in 1..n
   //---      x[i,j] in {0,1}, i in 1..n, j in 1..l[i]
   //---
   const MMKP_Param & appParam = m_model->getParam();
   UtilPrintFuncBegin(m_osLog, m_classTag,
		      "APPisUserFeasible()", appParam.LogLevel, 2);


   int            c, i, j, k, ij;
   bool           isFeasible = true;
   const MMKP_Instance & instance = m_model->getInstance();
   int            nGroupRows      = instance.getNGroupRows();
   int            nKnapRows       = instance.getNKnapRows();
   const double * capacity        = instance.getCapacity();
   const double * const * weight  = instance.getWeight();
   vector<int>    numInGroup(nGroupRows, 0);
   vector<double> weightInKnap(nKnapRows, 0.0);

   for(c = 0; c < n_cols; c++){
      CoinAssertDebug(UtilIsIntegral(x[c], tolZero));
      CoinAssertDebug(x[c] > (0.0 - tolZero));
      CoinAssertDebug(x[c] < (1.0 + tolZero));
      if(x[c] > 0.5){
	 pair<int,int> p = instance.getIndexInv(c);
	 i  = p.first;
	 j  = p.second;
	 ij = instance.getIndexIJ(i,j);
	 numInGroup[i]++;
	 for(k = 0; k < nKnapRows; k++){
	    weightInKnap[k] += weight[k][ij];
	 }
      }
   }
   for(i = 0; i < nGroupRows; i++){
      //printf("APPisUserFeasible numInGroup[%d]: %d", i, numInGroup[i]);
      if(numInGroup[i] != 1){
	 //printf(" --> NOT FEASIBLE");
	 isFeasible = false;
      }
      //printf("\n");
   }
   for(k = 0; k < nKnapRows; k++){
      //printf("APPisUserFeasible weightInKnap[%d]: %g, cap: %g",
      //     k, weightInKnap[k], capacity[k]);
      if(weightInKnap[k] >= (capacity[k] + tolZero)){
	 //printf(" --> NOT FEASIBLE");
	 isFeasible = false;
      }
      //printf("\n");
   }   

   //printf("APPisUserFeasible = %d\n", isFeasible);
   UtilPrintFuncEnd(m_osLog, m_classTag,
                    "APPisUserFeasible()", appParam.LogLevel, 2);

   return isFeasible;
}
//#endif

