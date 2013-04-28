//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2013, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "UtilMacros.h"
#include "ATM_Instance.h"


//===========================================================================//
void ATM_Instance::readInstance(string & fileNameA,
				string & fileNameD,
				string & fileNameAD){

   string   atm, date;
   int      a, d, ad, n_ad;
   double   K, B;
   
   ifstream isA, isD, isAD;
   ifstream isA2, isD2; //MSVS bug?
   char      dummy[10000];
   const int maxLine = 10000;

   //---
   //--- File format (.../Decomp/data/ATM)
   //---
   //--- dataA.txt:
   //---    a K[a]
   //--- dataD.txt:
   //---    d B[d]
   //--- dataAD.txt:
   //---    a d a[a,d] b[a,d] c[a,d] d[a,d] e[a,d]
   //---
   
   //---
   //--- open file streams
   //---
   UtilOpenFile(isA, fileNameA.c_str());
   UtilOpenFile(isD, fileNameD.c_str());
   UtilOpenFile(isAD, fileNameAD.c_str());

   //---
   //--- get number of atms
   //---
   printf("Reading %s\n", fileNameA.c_str());
   m_nAtms = 0;
   isA.getline(dummy, maxLine);
   while(!isA.eof()){
      isA >> atm;
      if(isA.eof())
	 break;
      isA >> K;
      m_nAtms++;
   }
   isA.close();

   //---
   //--- get number of dates
   //---
   printf("Reading %s\n", fileNameD.c_str());
   m_nDates = 0;
   isD.getline(dummy, maxLine);
   while(!isD.eof()){
      isD >> date;
      if(isD.eof())
	 break;
      isD >> B;
      m_nDates++;
   }
   isD.close();
   
   //---
   //--- allocate memory for storage of coefficients
   //---   open enough space as if dense
   //---
   n_ad   = m_nAtms * m_nDates;
   m_a_ad = new double[n_ad];
   m_b_ad = new double[n_ad];
   m_c_ad = new double[n_ad]; //(=-b)
   m_d_ad = new double[n_ad];
   m_e_ad = new double[n_ad];
   m_w_ad = new double[n_ad];
   m_B_d  = new double[m_nDates];
   m_K_a  = new double[m_nAtms];
   assert(m_a_ad && 
          m_b_ad && 
          m_c_ad && 
          m_d_ad && 
          m_e_ad &&
          m_w_ad &&
          m_B_d  && 
          m_K_a);
   
   //---
   //--- get data for atms
   //---
   UtilOpenFile(isA2, fileNameA.c_str());
   UtilOpenFile(isD2, fileNameD.c_str());
   m_nAtms = 0;
   isA2.getline(dummy, maxLine);
   while(!isA2.eof()){
      isA2 >> atm;
      if(isA2.eof())
	 break;
      m_strToIntAtms.insert(make_pair(atm, m_nAtms));
      m_intToStrAtms.push_back(atm);
      isA2 >> m_K_a[m_nAtms];
      m_nAtms++;
   }
   isA2.close();

   //---
   //--- get data for dates
   //---
   m_nDates = 0;
   isD2.getline(dummy, maxLine);
   while(!isD2.eof()){
      isD2 >> date;
      if(isD2.eof())
	 break;
      m_strToIntDates.insert(make_pair(date, m_nDates));
      m_intToStrDates.push_back(date);
      isD2 >> m_B_d[m_nDates];
      m_nDates++;
   }
   isD2.close();

   //---
   //--- get data for ATMS x DATES (we don't have data for all pairs)
   //---   
   printf("Reading %s\n", fileNameAD.c_str());
   map<string,int>::iterator mi;
   isAD.getline(dummy, maxLine);
   while(!isAD.eof()){
      isAD >> atm >> date;
      if(isAD.eof())
	 break;

      //get a,d index for this pair
      mi = m_strToIntAtms.find(atm);
      if(mi == m_strToIntAtms.end()){
	 printf("ERROR atm not found: %s\n", atm.c_str());
      }
      assert(mi != m_strToIntAtms.end());
      a = mi->second;

      mi = m_strToIntDates.find(date);
      if(mi == m_strToIntDates.end()){
	 printf("ERROR dates not found: %s\n", date.c_str());
      }
      assert(mi != m_strToIntDates.end());
      d = mi->second;

      ad = getIndexAD(a,d);
      m_pairsAD.push_back(ad);

      isAD >> m_a_ad[ad];
      isAD >> m_b_ad[ad];
      isAD >> m_c_ad[ad];
      isAD >> m_d_ad[ad];
      isAD >> m_e_ad[ad];
      isAD >> m_w_ad[ad];

      //printf("ad=%d atm=%s date=%s a=%g b=%g c=%g d=%g e=%g w=%g\n",
      //	     ad, atm.c_str(), date.c_str(), 
      //     m_a_ad[ad],
      //     m_b_ad[ad],
      //     m_c_ad[ad],
      //     m_d_ad[ad],
      //     m_e_ad[ad],
      //     m_w_ad[ad]);
   }
   isAD.close();

   printf("Number of ATMS  = %d\n", getNAtms());
   printf("Number of Dates = %d\n", getNDates());
   printf("Number of Pairs = %d\n", getNPairs());
   
}

//===========================================================================//
void ATM_Instance::generateRandom(const int nAtms,
				  const int nDates,
				  const int seed){

   /*
     Data from original:
        \\ordsrv3\ormp\sas\ATM_Badshah\atm_20ATMS_3\atm_doc
	nDates=272, nAtms=20, nPairs=4730 (max=5440)

     proc means data=FTPLIB.amul_atms_dates;
     var withdrawal allocation NET_IMPACT_AVG 
     NET_IMPACT_STD NORMAL_AVG NORMAL_STD TS1 TS2;
     run;
     
     The MEANS Procedure     
     Variable          Label                  Std Dev            Mean
     ................................................................
     WITHDRAWAL        WITHDRAWAL          1456368.37      1457077.41
     ALLOCATION        ALLOCATION          1752334.72      2068196.66
     NET_IMPACT_AVG    NET_IMPACT_AVG       0.8990607       1.1961954
     NET_IMPACT_STD    NET_IMPACT_STD       1.8979644       1.4240460
     NORMAL_AVG        NORMAL_AVG          1352731.38      1440849.71
     NORMAL_STD        NORMAL_STD           352658.50       364123.38
     TS1               TS1                 1267244.80      1371637.24
     S2                TS2                 1246864.33      1361954.95
     ................................................................
          

     Variable          Label                  Minimum         Maximum
     ................................................................
     WITHDRAWAL        WITHDRAWAL             8000.00      7080400.00
     ALLOCATION        ALLOCATION           100000.00      7020000.00
     NET_IMPACT_AVG    NET_IMPACT_AVG       0.0053384      18.7119586
     NET_IMPACT_STD    NET_IMPACT_STD       0.0046809      54.0086478
     NORMAL_AVG        NORMAL_AVG            38864.52      4375539.71
     NORMAL_STD        NORMAL_STD            26833.85      1141006.06
     TS1               TS1                   25245.45      5250885.71
     TS2               TS2                700.0000000      4182207.14
     ................................................................
     
     for{<a,d> in ATMS_DATES_THIS} do;
     CA[a,d] = (normal_avg[a,d]  * net_impact_avg[a,d] - ts_period_2[a,d]); 
     CB[a,d] = (ts_period_1[a,d] - ts_period_2[a,d]);
     CC[a,d] = (ts_period_2[a,d] - ts_period_1[a,d]);
     CD[a,d] = (normal_std[a,d]  * net_impact_std[a,d]); 
     CE[a,d] = (-actual_withdrawal[a,d] + ts_period_2[a,d]);
     end;

     These numbers are annoying big and causing lots of numerical 
     round-off issues. So, Let's scale by 1000.
   */
#define STDD 0
#define MEAN 1
#define MIN  2
#define MAX  3
   double s_withdrawal[4]   = {1456368,   1457077,   8000,       7080400};
   double s_allocation[4]   = {1752334,   2068196,   100000,     7020000};
   double s_netimpactAve[4] = {0.8990607, 1.1961954, 0.0053384, 18.7119586};
   double s_netimpactStd[4] = {1.8979644, 1.4240460, 0.0046809, 54.0086478};
   double s_normalAve[4]    = {13527318,  1440849,   38864,     4375539.71};
   double s_normalStd[4]    = {352658,    364123,    26833,     1141006};
   double s_ts1[4]          = {1267244,   1371637,   25245,     5250885};
   double s_ts2[4]          = {1246864,   1361954,   700,       4182207};
   double scale             = 1000;
   int i;
   for(i = 0; i < 4; i++){
      s_withdrawal[i] /= scale;
      s_allocation[i] /= scale;
      s_normalAve[i]  /= scale;
      s_normalStd[i]  /= scale;
      s_ts1[i]        /= scale;
      s_ts2[i]        /= scale;
   }


   int      nAD          = nAtms * nDates;
   double * withdrawal   = new double[nAD];
   double * allocation   = new double[nAD];
   double * netimpactAve = new double[nAD];
   double * netimpactStd = new double[nAD];
   double * normalAve    = new double[nAD];
   double * normalStd    = new double[nAD];
   double * ts1          = new double[nAD];
   double * ts2          = new double[nAD];
   assert(withdrawal   && allocation   && netimpactAve &&
	  netimpactStd && normalAve    && normalStd    &&
	  ts1          && ts2);
   
   string   fileNameAD   = "atm_randAD_";
   fileNameAD += UtilIntToStr(nAtms) + "_";
   fileNameAD += UtilIntToStr(nDates) + "_";
   fileNameAD += UtilIntToStr(seed) + ".txt";
   string   fileNameA   = "atm_randA_";
   fileNameA += UtilIntToStr(nAtms) + "_";
   fileNameA += UtilIntToStr(nDates) + "_";
   fileNameA += UtilIntToStr(seed) + ".txt";
   string   fileNameD   = "atm_randD_";
   fileNameD += UtilIntToStr(nAtms) + "_";
   fileNameD += UtilIntToStr(nDates) + "_";
   fileNameD += UtilIntToStr(seed) + ".txt";

   ofstream osAD, osA, osD;
   UtilOpenFile(osAD, fileNameAD.c_str());
   UtilOpenFile(osA,  fileNameA.c_str());
   UtilOpenFile(osD,  fileNameD.c_str());
   
   int a, d;
   srand(seed);
   //---
   //--- generate 'raw data' in N[mean,std-dev]
   //---
   int index = 0;//a * nDates + d
   for(a = 0; a < nAtms; a++){
      for(d = 0; d < nDates; d++){
	 do{
	    withdrawal[index]   = UtilNormRand(s_withdrawal[MEAN]  , 
					       s_withdrawal[STDD]);
	 }while( withdrawal[index] < s_withdrawal[MIN] ||
		 withdrawal[index] > s_withdrawal[MAX]);
	 do{
	 allocation[index]   = UtilNormRand(s_allocation[MEAN]  , 
					    s_allocation[STDD]);
	 }while( allocation[index] < s_allocation[MIN] ||
		 allocation[index] > s_allocation[MAX]);
	 do{
	    netimpactAve[index] = UtilNormRand(s_netimpactAve[MEAN], 
					    s_netimpactAve[STDD]);
	 }while( netimpactAve[index] < s_netimpactAve[MIN] ||
		 netimpactAve[index] > s_netimpactAve[MAX]);
	 do{
	    netimpactStd[index] = UtilNormRand(s_netimpactStd[MEAN], 
					       s_netimpactStd[STDD]);
	 }while( netimpactStd[index] < s_netimpactStd[MIN] ||
		 netimpactStd[index] > s_netimpactStd[MAX]);
	 do{
	    normalAve[index]    = UtilNormRand(s_normalAve[MEAN]   , 
					       s_normalAve[STDD]);
	 }while( normalAve[index] < s_normalAve[MIN] ||
		 normalAve[index] > s_normalAve[MAX]);
	 do{
	    normalStd[index]    = UtilNormRand(s_normalStd[MEAN]   , 
						s_normalStd[STDD]);
	 }while( normalStd[index] < s_normalStd[MIN] ||
		 normalStd[index] > s_normalStd[MAX]);
	 do{
	    ts1[index]          = UtilNormRand(s_ts1[MEAN]         , 
					       s_ts1[STDD]);
	 }while( ts1[index] < s_ts1[MIN] ||
		 ts1[index] > s_ts1[MAX]);
	 do{
	    ts2[index]          = UtilNormRand(s_ts2[MEAN]         , 
					    s_ts2[STDD]);
	 }while ( ts2[index] < s_ts2[MIN] ||
		  ts2[index] > s_ts2[MAX]);
	 index++; 
      }
   }
		
   //---
   //--- generate coefficients
   //--- 
   //CA[a,d] = (normal_avg[a,d]  * net_impact_avg[a,d] - ts_period_2[a,d]); 
   //CB[a,d] = (ts_period_1[a,d] - ts_period_2[a,d]);
   //CC[a,d] = (ts_period_2[a,d] - ts_period_1[a,d]);
   //CD[a,d] = (normal_std[a,d]  * net_impact_std[a,d]); 
   //CE[a,d] = (-actual_withdrawal[a,d] + ts_period_2[a,d]);
   double * ca = new double[nAD];
   double * cb = new double[nAD];
   double * cc = new double[nAD];
   double * cd = new double[nAD];
   double * ce = new double[nAD];
   assert(ca && cb && cc && cd && ce);
   
   index = 0;
   osAD << "a\td\tCA\tCB\tCC\tCD\tCD\tCE\tCW\n";
   for(a = 0; a < nAtms; a++){
      for(d = 0; d < nDates; d++){
	 ca[index] = normalAve[index] * netimpactAve[index] - ts2[index];
	 cb[index] = ts1[index] - ts2[index];
	 cc[index] = -cb[index];
	 cd[index] = normalStd[index] * netimpactStd[index];
	 ce[index] = -withdrawal[index] + ts2[index];
	 osAD << "ATM"  << UtilIntToStr(a) << "\t"
	      << "DATE" << UtilIntToStr(d) << "\t"
	      << setw(10) << UtilDblToStr(ca[index],0)
	      << setw(10) << UtilDblToStr(cb[index],0)
	      << setw(10) << UtilDblToStr(cc[index],0)
	      << setw(10) << UtilDblToStr(cd[index],0)
	      << setw(10) << UtilDblToStr(ce[index],0) 
	      << setw(10) << UtilDblToStr(withdrawal[index],0) 
	      << "\n";
	 index++;
      }
   }

   //---
   //--- generate B and K
   //---
   //--- f(a,d) = ca*x1[a] + cb*x2[a] + cc*x1[a]*x2[a] + cd*x3[a] + ce
   //---      x1,x2 in {0,1}, x3 >= 0
   //---
   //--- sum{a} f(a,d)           <= B[d], for d
   //--- |{d in D | f(a,d) <= 0} <= K[a], for a
   //---
   double x01, x1, x2, x3;
   double * f = new double[nAD];
   assert(f);
   index = 0;
   for(a = 0; a < nAtms; a++){
      x01 = UtilURand(0.0,1.0);
      x1  = x01 >= 0.5 ? 1 : 0;
      x01 = UtilURand(0.0,1.0);
      x2  = x01 >= 0.5 ? 1 : 0;
      x3 = UtilURand(0.0,1.0);
      printf("x1=%g x2=%g x3=%g\n", x1, x2, x3);
      for(d = 0; d < nDates; d++){
	 f[index]  = ca[index] * x1;
	 f[index] += cb[index] * x2;
	 f[index] += cc[index] * x1 * x2;
	 f[index] += cd[index] * x3;
	 f[index] += ce[index];	    
	 index++;
      }
   }   
   
   double * B    = new double[nDates];
   double   maxB = -1e20;
   for(d = 0; d < nDates; d++){
      B[d] = 0;
      for(a = 0; a < nAtms; a++){
	 B[d] += f[a * nDates + d];
      }
      if(B[d] > maxB) maxB=B[d];
   }
   //---
   //--- B=budget for cash flow
   //---   if negative does not make sense
   //---   protect against this
   //---
   osD << "d\tB\n";
   for(d = 0; d < nDates; d++){
      if(B[d] < 0)
	 B[d] = maxB / 2.0;
      osD << "DATE" << UtilIntToStr(d) << "\t"
	  << setw(10) << UtilDblToStr(B[d],0) << endl;
   }
 
   int * K    = new int[nAtms];
   int   maxK = 0;
   index = 0;
   for(a = 0; a < nAtms; a++){
      K[a] = 0;
      for(d = 0; d < nDates; d++){
	 K[a] += f[index] <= 0 ? 1 : 0;
	 index++;
      }
      if(K[a] > maxK) maxK=K[a];
   }
   //---
   //--- randomize it (and tighten constraint)
   //---
   osA << "a\tK\n";
   for(a = 0; a < nAtms; a++){
      //K[a] -= UtilURand(1, maxK/4);
      osA << "ATM" << UtilIntToStr(a) << "\t"
	  << setw(10) << K[a] << endl;
   }
   
   osAD.close();
   osA.close();
   osD.close();
   
   UTIL_DELARR(withdrawal);
   UTIL_DELARR(allocation);
   UTIL_DELARR(netimpactAve);
   UTIL_DELARR(netimpactStd);
   UTIL_DELARR(normalAve);
   UTIL_DELARR(normalStd);
   UTIL_DELARR(ts1);
   UTIL_DELARR(ts2);
   UTIL_DELARR(ca);
   UTIL_DELARR(cb);
   UTIL_DELARR(cc);
   UTIL_DELARR(cd);
   UTIL_DELARR(ce);   
   UTIL_DELARR(f);
   UTIL_DELARR(B);
   UTIL_DELARR(K);
}

