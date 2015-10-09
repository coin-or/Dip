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
// Copyright (C) 2002-2015, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//

#ifndef UTIL_PARAMETERS_INCLUDED
#define UTIL_PARAMETERS_INCLUDED

//===========================================================================//
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

//===========================================================================//
class UtilParameters {
private:
   std::map<std::string, std::string> m_paramMap;

public:
   UtilParameters()
      : m_paramMap() {}

   UtilParameters(int&   argc,
                  char* argv[]) :
      m_paramMap() {
      ScanCmdLineArgs(argc, argv);
   };

   ~UtilParameters() {};

   const std::map<std::string, std::string> getParamMap() { return m_paramMap; }
      
   void ScanCmdLineArgs(int&   argc,
                        char* argv[]);
   void   LoadParamFile(std::string& paramFileName);
   void   Add(std::string& section,
              std::string& name,
              std::string& value);
   void   Add(const char* section,
              const char* name,
              const char* value);
   std::string GetSetting(const char* name,
                          const char* defaultValue,
                          const char* section = NULL);
   int    GetSetting(const char* name,
                     const int    defaultValue,
                     const char* section = NULL);
   bool   GetSetting(const char* name,
                     const bool   defaultValue,
                     const char* section = NULL);
   long   GetSetting(const char* name,
                     const long   defaultValue,
                     const char* section = NULL);
   double GetSetting(const char* name,
                     const double defaultValue,
                     const char* section = NULL);

   std::string GetSetting(const char* name,
                          const std::string defaultValue,
                          const char* section = NULL) ;
   
   std::vector<std::string> &split(const std::string &s,
				   std::vector<std::string> &elems,
				   char delim = '@') {
      std::stringstream ss(s);
      std::string item;
      while (std::getline(ss, item, delim)) {
	 elems.push_back(item);
      }
      return elems;
   }

 private:
   
   std::string*     Find(const char* section,
                         const char* name);
};

#endif
