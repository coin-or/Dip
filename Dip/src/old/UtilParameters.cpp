//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//

//===========================================================================//
#include "UtilMacros.h"
#include "UtilParameters.h"

#include <cassert>
#include <iostream>
using namespace std;

#define MAXLINE 1024
//===========================================================================//

//TODO: need a template for usage, need application to register this, defaults
//TODO: this is an ugly combination of old C code and C++, C++-ify it

// ------------------------------------------------------------------------- //
void UtilParameters::ScanCmdLineArgs(int&   argc,
                                     char* argv[])
{
   int i, j;

   //---
   //--- if there are no arguments, return
   //---
   if (argc == 0 || argv == NULL) {
      return;
   }

   j = 1;
   string paramFileName;

   for (i = 1; i < argc; i++) {
      //---
      //--- if "--param" is specified
      //---
      if (strcmp(argv[i], "--param") == 0) {
         //---
         //---   if a filename follows the flag
         //---           grab the filename
         //---
         if (((i + 1) < argc)
               && (argv[i+1][0] != '-' || argv[i+1][1] != '-')) {
            paramFileName = argv[++i];
         }

         continue;
      }

      argv[j++] = argv[i];
   }

   argc = j;
   //---
   //--- load the parameter file
   //---
   LoadParamFile(paramFileName);
   //---
   //--- enter the command line flags into the parameter map
   //---   format is --SECTION:PARAMETER value
   //---
   char cmdBuf[MAXLINE];

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-' && argv[i][1] == '-') {
         string name (argv[i] + 2);
         string value("");
         strcpy(cmdBuf, argv[i] + 2);
         char* ptr1    = strtok(cmdBuf, ":");
         char* ptr2    = strtok(NULL, ":");
         char* section = NULL;
         char* parm    = NULL;
         printf("\nptr1 = %s", ptr1);
         printf("\nptr2 = %s", ptr2);
         fflush(stdout);

         if (ptr2 == NULL) {
            //---
            //--- section is NULL
            //---
            parm = ptr1;
         } else {
            //---
            //--- section is not NULL
            //---
            section = ptr1;
            parm    = ptr2;
         }

         if (((i + 1) < argc) && (argv[i+1][0] != '-' || argv[i+1][1] != '-')) {
            value = argv[++i];
         }

         Add(section, parm, value.c_str());
         continue;
      }
   }

   //STOP - what is the point of this?
   //
   //      -----   remove all the flags from the command line
   //
   j = 1;

   for (i = 1; i < argc; i++) {
      if (argv[i][0] == '-' && argv[i][1] == '-') {
         string name (argv[i] + 2);
         string value("");

         if (((i + 1) < argc)
               && (argv[i+1][0] != '-' || argv[i+1][1] != '-')) {
            value = argv[++i];
         }

         continue;
      }

      argv[j++] = argv[i];
   }

   argc = j;
}

// ------------------------------------------------------------------------- //
void UtilParameters::LoadParamFile(string& paramFileName)
{
   char     buf[MAXLINE];
   char*    ptr = NULL;
   string   curSection("");
   string   bufStr    ("");
   string   name      ("");
   string   value     ("");
   //---
   //--- open the stream pointer
   //---
   ifstream is(paramFileName.c_str());

   if (!is) {
      return;
   }

   //---
   //---   foreach line in the file
   //---     skip comments (#) and blank lines
   //---
   int lineNum = 0;

   while (!is.eof()) {
      is.getline(buf, sizeof(buf) - 1);

      if (is.eof()) {
         break;
      }

      lineNum++;
      ptr = strchr(buf, '#');

      if (ptr != NULL) {
         *ptr = '\0';
      }

      //TODO: move all to use string? do we need buf?
      bufStr = buf;
      bufStr = UtilStrTrim(bufStr);
      strcpy(buf, bufStr.c_str());

      if (strlen(buf) < 3) {
         continue;
      }

      //---
      //--- if line is '[section]'
      //---   create a new section
      //---
      if (buf[0] == '[') {
         ptr = strchr(buf + 1, ']');

         if (ptr == NULL) {
            cerr << "UtilParameters: syntax error on line "
                 << lineNum << " '" << buf << "'" << endl;
         }

         *ptr = '\0';
         curSection = buf + 1;
         continue;
      }

      //---
      //--- if line is 'name = value'
      //---   create a new name/value pair in the current section
      //---
      ptr = strchr(buf, '=');

      if (ptr != NULL) {
         *ptr++ = '\0';
      }

      name  = buf;
      value = "1";
      name  = UtilStrTrim(name);
      name  = UtilStrToLower(name);

      if (ptr != NULL) {
         //TODO: look into old code -> value=Expand(ptr)
         value = ptr;
         value = UtilStrTrim(value);
      }

      Add(curSection, name, value);
   }

   //---
   //--- close file stream
   //---
   is.close();
}

// ------------------------------------------------------------------------- //
void UtilParameters::Add(string& sSection,
                         string& sName,
                         string& sValue)
{
   UtilParam utilParam;
   string    keyname ("");
   keyname  = UtilStrToLower(UtilStrTrim(sSection));
   keyname += "@";
   keyname += UtilStrToLower(UtilStrTrim(sName));
   utilParam.paramName = UtilStrTrim(sValue);
   utilParam.isUsed    = false;
   //TODO: why doesn't insert override??
   // m_paramMap.insert(make_pair(keyname, utilParam));
   m_paramMap[keyname].paramName = utilParam.paramName;
   m_paramMap[keyname].isUsed    = utilParam.isUsed;
}

// ------------------------------------------------------------------------- //
void UtilParameters::Add(const char* section,
                         const char* name,
                         const char* value)
{
   UtilParam utilParam;
   string    keyname ("");
   string    sSection("");
   string    sName   (name);
   string    sValue  (value);

   if (section) {
      sSection = section;
      keyname  = UtilStrToLower(UtilStrTrim(sSection));
   }

   keyname += "@";
   keyname += UtilStrToLower(UtilStrTrim(sName));
   utilParam.paramName = UtilStrTrim(sValue);
   utilParam.isUsed    = false;
   //TODO: why doesn't insert override??
   // m_paramMap.insert(make_pair(keyname, utilParam));
   m_paramMap[keyname].paramName = utilParam.paramName;
   m_paramMap[keyname].isUsed    = utilParam.isUsed;
}

// ------------------------------------------------------------------------- //
UtilParam* UtilParameters::FindEntry(const char* section,
                                     const char* name)
{
   string    keyname ("");
   string    sSection("");
   string    sName   (name);

   if (section) {
      sSection = section;
      keyname  = UtilStrToLower(UtilStrTrim(sSection));
   }

   keyname += "@";
   keyname += UtilStrToLower(UtilStrTrim(sName));
   map<string, UtilParam>::iterator it;
   it = m_paramMap.find(keyname);

   if (it == m_paramMap.end()) {
      return NULL;
   } else {
      return &(it->second);
   }
}

// ------------------------------------------------------------------------- //
string* UtilParameters::Find(const char* section,
                             const char* name)
{
   UtilParam* utilParam = FindEntry(section, name);

   if (utilParam == NULL) {
      Add(section, name, "(undefined)");
      return NULL;
   } else {
      return &utilParam->paramName;
   }
}

// ------------------------------------------------------------------------- //
string UtilParameters::GetSetting(const char* name,
                                  const char* defaultValue,
                                  const char* section)
{
   //---
   //--- build the qualified name using the section
   //--- if the parameter is not found, return the default
   //--- else convert the string to the appropriate type
   //---
   string* pVal = Find(section, name);

   if (pVal == NULL) {
      return string(defaultValue);
   }

   return *pVal;
}

// ------------------------------------------------------------------------- //
int UtilParameters::GetSetting(const char* name,
                               const int    defaultValue,
                               const char* section)
{
   //---
   //--- build the qualified name using the section
   //--- if the parameter is not found, return the default
   //--- else convert the string to the appropriate type
   //---
   string* pVal = Find(section, name);

   if (pVal == NULL) {
      return defaultValue;
   }

   int value = atoi(pVal->c_str());
   return value;
}

// ------------------------------------------------------------------------- //
bool UtilParameters::GetSetting(const char* name,
                                const bool   defaultValue,
                                const char* section)
{
   //---
   //--- build the qualified name using the section
   //--- if the parameter is not found, return the default
   //--- else convert the string to the appropriate type
   //---
   string* pVal = Find(section, name);

   if (pVal == NULL) {
      return defaultValue;
   }

   bool value = atoi(pVal->c_str()) != 0;
   return value;
}

// ------------------------------------------------------------------------- //
long UtilParameters::GetSetting(const char* name,
                                const long   defaultValue,
                                const char* section)
{
   //---
   //--- build the qualified name using the section
   //--- if the parameter is not found, return the default
   //--- else convert the string to the appropriate type
   //---
   string* pVal = Find(section, name);

   if (pVal == NULL) {
      return defaultValue;
   }

   long value = atol(pVal->c_str());
   return value;
}

// ------------------------------------------------------------------------- //
double UtilParameters::GetSetting(const char* name,
                                  const double defaultValue,
                                  const char* section)
{
   //---
   //--- build the qualified name using the section
   //--- if the parameter is not found, return the default
   //--- else convert the string to the appropriate type
   //---
   string* pVal = Find(section, name);

   if (pVal == NULL) {
      return defaultValue;
   }

   char*    pEnd  = NULL;
   double   value = strtod(pVal->c_str(), &pEnd);
   return value;
}

#if 0
// ------------------------------------------------------------------------- //
void UtilParameters::Dump(std::ostream& os) const
{
   //
   //    -----   foreach section
   //
   SECTION_MAP::const_iterator iSection;

   for (iSection =  m_sectionMap.begin();
         iSection != m_sectionMap.end(); iSection++) {
      //
      //        -----           output the section name
      //
      std::string sectionName = iSection->first;

      if (sectionName.length() > 0) {
         os << std::endl;
         os << "[" << sectionName << "]" << std::endl;
      }

      const NV_MAP& nv = iSection->second;

      //
      //        -----           foreach parameter in the section
      //
      NV_MAP::const_iterator iParm;

      for (iParm = nv.begin(); iParm != nv.end(); iParm++) {
         const std::string& name = iParm->first;
         const std::string& value = iParm->second.second;
         //
         //    -----                   output the name/value pair
         //
         os << "[" << iParm->second.first << "] "
            << name << "    = "
            << value << std::endl;
      }
   }
}
#endif
