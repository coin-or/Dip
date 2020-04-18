#!/usr/bin/env python

from setuptools import setup, Extension, find_packages
import subprocess, os, sys
from os.path import join, dirname
from os import listdir

PROJECT = 'coinor.dippy'
VERSION = '1.95.2'
URL = 'https://github.com/coin-or/Dip'
DESC = u'DIP Python Interface'

def read_file(file_name):

    file_path = join(
        dirname(__file__),
        file_name
        )

    return open(file_path).read()

try:
    CoinDir = os.environ['COIN_INSTALL_DIR']
except:
    from os.path import abspath, dirname

    try:
        location = dirname(
            check_output(['which', 'clp']).strip()).decode('utf-8')
        CoinDir = abspath(join(location, ".."))
    except:
        pass
            
def get_libs():
    '''
    Return a list of distinct library names used by ``dependencies``.
    '''
    libs = []

    try:
        from subprocess import check_output
        
        flags = (check_output(['pkg-config', '--libs', 'dip'])
                 .strip().decode('utf-8'))
        libs = [flag[2:] for flag in flags.split()
                if flag.startswith('-l')]
        libDirs = [flag[2:] for flag in flags.split()
                   if flag.startswith('-L')]
        flags = (check_output(['pkg-config', '--cflags', 'dip'])
                 .strip().decode('utf-8'))
        incDirs = [flag[2:] for flag in flags.split() if
                   flag.startswith('-I')]

    except:
        if CoinDir != None:
            
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
                
        else:
            raise Exception('''
            Could not find location of COIN installation.
            Please ensure that either 
            * COIN_INSTALL_DIR is set to the location of the installation,
            * PKG_CONFIG_PATH points to the location of the .pc files, or
            * The cbc executable is in your executable path and is installed
            at the same location as the libraries. 
            ''')

    return libs, libDirs, incDirs

operatingSystem = sys.platform
if 'linux' in operatingSystem:
    operatingSystem = 'linux'
elif 'darwin' in operatingSystem:
    operatingSystem = 'mac'
elif 'win' in operatingSystem:
    operatingSystem = 'windows'

libs, libDirs, incDirs = get_libs()

macros = [('__DECOMP_LP_CLP__', None)]

files = ['DippyDecompAlgo.cpp',
         'DippyDecompApp.cpp',
         'DippySolve.cpp',
         'DippyModule.cpp',
         'DippyPythonUtils.cpp',
         ]

sources = [join('src', 'dippy', f) for f in files]

modules=[Extension('coinor.dippy._dippy', 
                   sources, 
                   libraries=libs,
                   include_dirs=incDirs,
                   library_dirs=libDirs,
                   define_macros=macros)]

cvpmp_instance_files = [join('Instances', f) for f in
                        listdir('src/dippy/examples/cvpmp/Instances')]

setup(name=PROJECT,
      version=VERSION,
      description=DESC,
      long_description=read_file('README.md'),
      long_description_content_type='text/markdown',
      author='''Michael O'Sullivan, Qi-Shan Lim, Stuart Mitchell''',
      author_email='michael.osullivan@auckland.ac.nz',
      maintainer='Ted Ralphs',
      maintainer_email='ted@lehigh.edu',
      url=URL,
      license='Eclipse Public License',
      namespace_packages=['coinor'],
      #There must be a better way
      packages=[pkg.replace('src','coinor') for pkg in find_packages()],
      package_dir = {'coinor': 'src'},
      package_data = {'coinor.dippy.examples.cvpmp': cvpmp_instance_files},
      install_requires=['pulp>=1.5.4','coinor.gimpy>=2.0.0'],
      ext_modules=modules
     )


