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

#ifndef UTIL_APP_INCLUDED
#define UTIL_APP_INCLUDED

#include <string>
#include <map>

#include "UtilParameters.h"

class UtilApp
{
public:
    UtilApp(int& argc, char* argv[]);
    UtilApp() 
      : m_parms(),
      m_machine(""),
      m_program(""),
      m_fullPathname(""),
      m_pid(0){};
      ~UtilApp();

    static UtilApp& TheApp() {return *m_theApp;}

    const std::string& Machine();
    const std::string& Program() const;
    const std::string& FullPathname() const;
    int         Pid();
  //float       RunTime();      // in seconds and fractions

    bool        GetSetting(const char* name,
                           bool defaultValue=true,
                           const char* section = NULL) const;
    std::string GetSetting(const char* name,
                           const char* defaultValue,
                           const char* section = NULL) const;
    short       GetSetting(const char* name,
                           short defaultValue,
                           const char* section = NULL) const;
    int         GetSetting(const char* name,
                           int defaultValue,
                           const char* section = NULL) const;
    long        GetSetting(const char* name,
                           long defaultValue,
                           const char* section = NULL) const;
/*
    float       GetSetting(const char* name,
                           float defaultValue,
                           const char* section = NULL) const;
*/
    double       GetSetting(const char* name,
                           double defaultValue,
                           const char* section = NULL) const;

    bool        GetSetting(const char* name,
                           bool defaultValue,
                           const std::string& section) const
        {return GetSetting(name,defaultValue,section.c_str());}
    
    std::string GetSetting(const char* name,
                           const char* defaultValue,
                           const std::string& section) const
        {return GetSetting(name,defaultValue,section.c_str());}
    short       GetSetting(const char* name,
                           short defaultValue,
                           const std::string& section) const
        {return GetSetting(name,defaultValue,section.c_str());}
    int         GetSetting(const char* name,
                           int defaultValue,
                           const std::string& section) const
        {return GetSetting(name,defaultValue,section.c_str());}
    long        GetSetting(const char* name,
                           long defaultValue,
                           const std::string& section) const
        {return GetSetting(name,defaultValue,section.c_str());}
    double       GetSetting(const char* name,
                           double defaultValue,
                           const std::string& section) const
        {return GetSetting(name,defaultValue,section.c_str());}

    UtilParameters m_parms;

 public:
    void        LoadParmFile(std::string& fname);

protected:
  //void        add_to_environ(const char* env_var);
private:
    UtilApp&     operator=(const UtilApp& rhs);
    bool        operator==(const UtilApp& rhs) const;
    bool        operator<(const UtilApp& rhs) const;
#if 0
    std::string QualifiedName(const std::string& name,
                              const std::string& section) const;
    std::string QualifiedName(const char* name, const char* section) const;
#endif
    //void        LoadParmFile(std::string& fname);
    void        ScanCmdLineArgs(int& argc, char* argv[]);
protected:
//	-----	process variables
    std::string m_machine;
    std::string m_program;
    std::string m_fullPathname;
    int         m_pid;
  //UtilTimer    m_runTime;

//	-----	parm file variables
    typedef std::map<std::string,std::string,less<std::string> > PARM_MAP;
//    PARM_MAP    m_parms;

    static UtilApp* m_theApp;
};


#endif
