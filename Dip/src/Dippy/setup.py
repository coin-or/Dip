#!/usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension

from os.path import join, dirname

PROJECT = 'dippy'
VERSION = '1.4.3'
URL = 'https://projects.coin-or.org/Dip/wiki/Dippy'
AUTHOR_EMAIL = u''
DESC = 'DIP Python Interface'


def _src(files):
    return [join('dippy', f) for f in files]

def read_file(file_name):
    file_path = join(
        dirname(__file__),
        file_name
        )
    return open(file_path).read()

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
      package_dir = {'.': '.'},
      install_requires=['coinor.pulp'],
      ext_modules=[Extension('_dippy', 
                             _src(['DippyDecompAlgo.cpp',
                              'DippyDecompApp.cpp', 
                              'DippySolve.cpp',
                              'DippyModule.cpp',
                              'DippyPythonUtils.cpp']),
                             libraries=['Decomp', 'Alps', 'CbcSolver', 'Cgl',
                                        'Cbc', 'OsiClp', 'Osi', 'OsiCbc', 'Clp',
                                        'CoinUtils', 'm'],
                             define_macros=[('__DECOMP_LP_CLP__', None)])],
     )


