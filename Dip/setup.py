#!/usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension
import subprocess, os, sys
from os.path import join, dirname

PROJECT = 'dippy'
VERSION = '1.4.3'
URL = 'https://projects.coin-or.org/Dip/wiki/Dippy'
AUTHOR_EMAIL = u''
DESC = 'DIP Python Interface'

coin_install_dir = os.environ['COIN_INSTALL_DIR']

def read_file(file_name):

    file_path = join(
        dirname(__file__),
        file_name
        )

    return open(file_path).read()

def get_libs(dir):
    '''
    Return a list of distinct library names used by ``dependencies``.
    '''
    with open(join(dir, 'share', 'coin',
                   'doc', 'Dip', 'dip_addlibs.txt')) as f:
        link_line = f.read()
        if operatingSystem == 'windows':
            libs = [flag[:-4] for flag in link_line.split() if
                    flag.endswith('.lib')]
        else:
            libs = [flag[2:] for flag in link_line.split() if
                    flag.startswith('-l')]
    return libs

operatingSystem = sys.platform
if 'linux' in operatingSystem:
    operatingSystem = 'linux'
elif 'darwin' in operatingSystem:
    operatingSystem = 'mac'
elif 'win' in operatingSystem:
    operatingSystem = 'windows'

try:
    coin_install_dir = os.environ['COIN_INSTALL_DIR']
except KeyError:
    raise Exception('Please set the environment variable COIN_INSTALL_DIR' +
                    'to the location of the COIN installation')

if len(coin_install_dir.split(';')) > 1:
    raise Exception('Error: More than one directory listed in COIN_INSTALL_DIR')

libraries = get_libs(coin_install_dir)

macros = [('__DECOMP_LP_CLP__', None)]

files = ['DippyDecompAlgo.cpp',
         'DippyDecompApp.cpp',
         'DippySolve.cpp',
         'DippyModule.cpp',
         'DippyPythonUtils.cpp',
         ]

sources = [join('src/dippy', f) for f in files]

lib_dirs = [join(coin_install_dir, 'lib', 'intel'),
            join(coin_install_dir, 'lib')]

modules=[Extension('_dippy', 
                   sources, 
                   libraries=libraries,
                   include_dirs=[join(coin_install_dir, 'include', 'coin')],
                   library_dirs=lib_dirs,
                   define_macros=macros)]

setup(name=PROJECT,
      version=VERSION,
      description=DESC,
      long_description=read_file('README.rst'),
      author=read_file('AUTHORS'),
      author_email=AUTHOR_EMAIL,
      url=URL,
      license=read_file('LICENSE'),
#      namespace_packages=['coinor'],
      packages=['dippy'],
      package_dir = {'': 'src'},
      install_requires=['coinor.pulp'],
      ext_modules=modules
     )


