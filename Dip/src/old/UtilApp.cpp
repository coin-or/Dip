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

#ifndef WIN32
#include <unistd.h>     //linux
#include <sys/utsname.h>
#include <libgen.h>     //linux
#endif


//#include <stdio.h>      //linux
//#include <stdlib.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <cstring>
using namespace std;

#include "UtilApp.h"


 UtilApp* UtilApp::m_theApp = NULL;




UtilApp::UtilApp(int& argc, char* argv[])
: m_parms(),
  m_machine(""),
  m_program(""),
  m_fullPathname(""),
  m_pid(0){
  m_theApp = this;
//
//      -----   initialize the base app
//
    ScanCmdLineArgs(argc,argv);
}

UtilApp::~UtilApp()
{
}

void UtilApp::ScanCmdLineArgs(int& argc, char* argv[])
{
    int ii, jj;

//
//----- if there are no arguments
//-----         return
//
    if(argc == 0 || argv == NULL)
        return;
//
//      -----   grab the executable name from argv[0] and strip off
//      -----   the ".exe" suffix (if it is there)
//
    m_fullPathname = argv[0];
#ifdef CYGWIN
  m_program = "cygwinApp";
#endif
#ifdef WIN32
  m_program = "winApp";
#endif
#ifndef CYGWIN
#ifndef WIN32
  m_program = basename(const_cast<char*>(m_fullPathname.c_str()));
#endif
#endif
    int len = static_cast<int>(m_program.length());
    if(len > 4)
    {
        char *suffix = const_cast<char*>(m_program.c_str())+len-4;
        if(strcmp(suffix,".exe")==0)
        {
            *suffix = '\0';
            std::string tmpS = m_program.c_str();
            m_program = tmpS;
        }
    }

    jj = 1;

    std::string parmFileName;
    for(ii=1;ii<argc;ii++)
    {
//
//      -----   if "-dir" is specified
//      -----           change the current directory
//
        if(strcmp(argv[ii],"--dir")==0)
        {
            if(ii+1 < argc && (argv[ii+1][0] != '-' || argv[ii+1][1] != '-'))
            {
                ii++;
                chdir(argv[ii]);//work in linux?
            }
            continue;
        }
//
//      -----   if "-parm" is specified
//
        if(strcmp(argv[ii],"--parm")==0)
        {
//
//      -----   if a filename follows the flag
//      -----           grab the filename
//      -----           also set the MBG_PARM environment variable
//      ----    else default the name based on the program name
//
            if(ii+1 < argc && (argv[ii+1][0] != '-' || argv[ii+1][1] != '-')) 
            {
                parmFileName = argv[++ii];
                char buf[1024];
                sprintf(buf,"MBG_PARM=%s",parmFileName.c_str());
                //add_to_environ(buf);
            }
            continue;
        }
        argv[jj++] = argv[ii];
    }
    argc = jj;
//
//      -----   if the parameter file was not specified
//      -----           default to the name of the program + .parm
//
    if(parmFileName.empty())
        parmFileName = m_program + ".parm";

    LoadParmFile(parmFileName);
//
//      -----   enter the command line flags into the parameter map
//
    for(ii=1;ii<argc;ii++)
    {
        if(argv[ii][0] == '-' && argv[ii][1] == '-')
        {
          char * next_token;
          std::string name = argv[ii]+2;
          std::string value;
#if 1
          char cmdBuf[100];
          strcpy(cmdBuf,argv[ii]+2);
          //#ifdef WIN32
          //char* ptr1 = strtok_s(cmdBuf,":",&next_token);
          //char* ptr2 = strtok_s(NULL,":",&next_token);
          //#else          
          char * ptr1 = strtok(cmdBuf,":");
          char * ptr2 = strtok(NULL,":");
          //#endif
          char* section = NULL;
          char* parm;
          if(ptr2 == NULL)
            parm = ptr1;
          else 
            {
                section = ptr1;
                parm = ptr2;
            }
            if(ii+1 < argc && (argv[ii+1][0] != '-' || argv[ii+1][1] != '-'))
                value = argv[++ii];
            else 
                value = "1";

            m_parms.Add(section,parm,value.c_str());
#endif
            continue;
        }
    }
//
//      -----   remove all the flags from the command line
//
    jj = 1;
    for(ii=1;ii<argc;ii++)
    {
        if(argv[ii][0] == '-' && argv[ii][1] == '-')
        {
            std::string name = argv[ii]+2;
            std::string value;

            if(ii+1 < argc && (argv[ii+1][0] != '-' || argv[ii+1][1] != '-'))
            {
                value = argv[++ii];
            }
            continue;
        }
        argv[jj++] = argv[ii];
    }
    argc = jj;
}

void UtilApp::LoadParmFile(std::string& fname)
{
//
//      -----   open the parm file
//
    std::ifstream ifParm(fname.c_str());
    if(!ifParm)
        return;
#if 1
    m_parms.Load(ifParm);
#else
//
//      -----   read it
//
    std::string     curSection;
//
//      -----   foreach line in the file
//      -----           skip comments and blank lines
//
    char buf[1024];
    int lineNum = 0;
    while(!ifParm.eof())
    {
        ifParm.getline(buf,sizeof(buf)-1);
        if(ifParm.eof())
            break;

        lineNum++;
        trim(buf);
        if(strlen(buf) < 3)
            continue;
//
//      -----           if line is '[section]'
//      -----                   create a new section
//
        if(buf[0] == '[')
        {
            char* ptr = strchr(buf+1,']');
            if(ptr == NULL)
            {
                std::cerr << fname << ": syntax error on line " << lineNum
                          << " '" << buf << "'" << endl;
            }
            *ptr = '\0';
            curSection = buf+1;
            continue;
        }
//
//      -----           if line is 'name = value'
//      -----                   create a new name/value pair in the
//      -----                   current section
//
        bool bEqual = strchr(buf,'=') != NULL;

        const char* ptr = strtok(buf,"[ \t=]");
        std::string name(ptr);
        trim(name);
        lower(name);

        std::string parm = QualifiedName(name,curSection);

        if(bEqual)
        {
            ptr = strtok(NULL,"=");
            if(ptr == NULL)
                continue;
            ptr = strtok(NULL,"[=#;\n]");
            if(ptr == NULL)
                continue;
            std::string value(ptr);
            trim(value);

            m_parms[parm] = value;

            continue;
        } else { // no "=", must be a bool, so it is true just by being there
            m_parms[parm] = "1";
        }
//
//      -----           all other lines are errors
//
        cerr << fname << ": invalid parm line '" << buf << "'" << endl;
    }
#endif
    ifParm.close();
}

bool UtilApp::GetSetting(const char* name, bool unUsed, const char* section) const
{
//
//      -----   build the qualified name using the section
//      -----   if the parameter is not found, return the default
//      -----   else convert the string to the appropriate type
//
    return m_parms.GetSetting(name,unUsed,section);
}

std::string UtilApp::GetSetting(const char* name,
                               const char* defaultValue,
                               const char* section) const
{
//
//      -----   build the qualified name using the section
//      -----   if the parameter is not found, return the default
//      -----   else convert the string to the appropriate type
//
    return m_parms.GetSetting(name,defaultValue,section);
}

short UtilApp::GetSetting(const char* name,
                         short defaultValue,
                         const char* section) const
{
//
//      -----   build the qualified name using the section
//      -----   if the parameter is not found, return the default
//      -----   else convert the string to the appropriate type
//
    return m_parms.GetSetting(name,defaultValue,section);
}

int UtilApp::GetSetting(const char* name,
                       int defaultValue,
                       const char* section) const
{
//
//      -----   build the qualified name using the section
//      -----   if the parameter is not found, return the default
//      -----   else convert the string to the appropriate type
//
    return m_parms.GetSetting(name,defaultValue,section);
}

long UtilApp::GetSetting(const char* name,
                        long defaultValue,
                        const char* section) const
{
//
//      -----   build the qualified name using the section
//      -----   if the parameter is not found, return the default
//      -----   else convert the string to the appropriate type
//
    return m_parms.GetSetting(name,defaultValue,section);
}

double UtilApp::GetSetting(const char* name,
                         double defaultValue,
                         const char* section) const
{
//
//      -----   build the qualified name using the section
//      -----   if the parameter is not found, return the default
//      -----   else convert the string to the appropriate type
//
    return m_parms.GetSetting(name,defaultValue,section);
}

const std::string& UtilApp::Machine()
{
    if(!m_machine.empty())
        return m_machine;

#ifdef Linux
    m_machine = "PC";
#else

#ifndef UNIX 
    m_machine = "PC";
#else
    struct utsname info;
    extern int errno;

    if(uname(&info)<0)
    {
        char buf[100];
        sprintf(buf,"?%d?",errno);
        m_machine = buf;
    } else
        m_machine = info.nodename;
#endif
#endif

    return m_machine;
}

const std::string& UtilApp::Program() const
{
    return m_program;
}

const std::string& UtilApp::FullPathname() const
{
    return m_fullPathname;
}

int UtilApp::Pid()
{
    if(m_pid != 0)
        return m_pid;
#ifdef UNIX
    m_pid = getpid();
#else
    m_pid = 0;
#endif
    return m_pid;
}
