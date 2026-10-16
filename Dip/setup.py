#!/usr/bin/env python

from setuptools import setup, Extension, find_namespace_packages
import subprocess, os, sys
from os.path import join, dirname
from os import listdir
from subprocess import check_output

def read_file(file_name):

    file_path = join(
        dirname(__file__),
        file_name
        )

    return open(file_path).read()

CoinDir = None

libs = []
libDirs = []
incDirs = []

try:
    if len(sys.argv) > 1 and (sys.argv[1] == "sdist" or sys.argv[1] == "egg_info"):
        # Do not need CoinDir
        pass
    else:
        CoinDir = os.environ['COIN_INSTALL_DIR']
except:
    # If user didn't supply location, then try pkg-config
    try:
        for p in ['dip','alps','cbc','cgl','osi-clp','clp','osi','coinutils']:
            flags = (check_output(['pkg-config', '--libs', p])
                     .strip().decode('utf-8'))
            for flag in flags.split():
                if flag.startswith('-l') and flag[2:] not in libs:
                    libs.append(flag[2:]) 
                if flag.startswith('-L') and flag[2:] not in libDirs:
                    libDirs.append(flag[2:]) 
            flags = (check_output(['pkg-config', '--cflags', p])
                     .strip().decode('utf-8'))
            for flag in flags.split():
                if flag.startswith('-I') and flag[2:] not in incDirs:
                    incDirs.append(flag[2:]) 
    except:
        #If pkg-config fails, then look for an installed Cbc
        try:
            location = dirname(
                check_output(['which', 'cbc']).strip()).decode('utf-8')
            CoinDir = abspath(join(location, ".."))
        except:
            #Otherwise, raise an exception
            raise Exception('''
            Could not find location of COIN installation.
            Please ensure that either 
            * COIN_INSTALL_DIR is set to the location of the installation,
            * PKG_CONFIG_PATH points to the location of the .pc files, or
            * The cbc executable is in your executable path and is installed
            at the same location as the libraries. 
            ''')

operatingSystem = sys.platform
if 'linux' in operatingSystem:
    operatingSystem = 'linux'
elif 'darwin' in operatingSystem:
    operatingSystem = 'mac'
elif 'win' in operatingSystem:
    operatingSystem = 'windows'

if CoinDir != None:
    # We come here if user supplied the installation directory or pkg-config failed
    if operatingSystem == 'windows':
        if os.path.exists(join(CoinDir, 'lib', 'Cbc.lib')):
            libs = ['Decomp', 'OsiSym', 'Sym', 'OsiCbc', 'CbcSolver', 'Cbc', 'Cgl',
                    'OsiClp', 'ClpSolver',  'Clp', 'Osi', 'Alps', 'CoinUtils']
        else:
            libs = ['libDecomp', 'libOsiSym', 'libSym', 'libOsiCbc', 'libCbcSolver',
                    'libCbc', 'libCgl', 'libOsiClp', 'libClpSolver',  'libClp', 'libOsi',
                    'libAlps', 'libCoinUtils']
    else:
        libs = ['Decomp', 'OsiSym', 'Sym', 'OsiCbc', 'CbcSolver', 'Cbc', 'Cgl',
                'OsiClp', 'ClpSolver',  'Clp', 'Osi', 'Alps', 'CoinUtils']
                
    libDirs = [join(CoinDir, 'lib')]
    incDirs = [join(CoinDir, 'include', 'coin')] 
            
incDirs.extend(join('src', 'coinor', 'dippy'))

macros = [('__DECOMP_LP_CLP__', None)]

files = ['DippyDecompAlgo.cpp',
         'DippyDecompApp.cpp',
         'DippySolve.cpp',
         'DippyModule.cpp',
         'DippyPythonUtils.cpp',
         ]

sources = [join('src', 'coinor', 'dippy', f) for f in files]

modules=[Extension('coinor.dippy._dippy', 
                   sources, 
                   libraries=libs,
                   include_dirs=incDirs,
                   library_dirs=libDirs,
                   define_macros=macros)]

cvpmp_instance_files = [join('Instances', f) for f in
                        listdir(join('src','coinor', 'dippy', 'examples',
                                     'cvpmp'))]

setup(packages=find_namespace_packages(where='src/', include=['coinor.dippy*']),
      package_dir = {'': 'src'},
      package_data = {'coinor.dippy.examples.cvpmp': cvpmp_instance_files},
      ext_modules=modules
     )


