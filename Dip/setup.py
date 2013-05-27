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

def get_libs():
       
    p = subprocess.Popen(['pkg-config', '--libs', 'dip'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout_output, stderr_output = p.communicate()

    if p.returncode != 0 :
        raise Exception('pkg-config terminated with status:'+
                        '%d. stderr follows: %s' 
                        % (p.returncode, stderr_output))
    elif stderr_output:
        print stderr_output

    lib_dirs = [flag[2:] for flag in stdout_output.split()
                if flag.startswith('-L')]

    libs =     [flag[2:] for flag in stdout_output.split()
                if flag.startswith('-l')]

    return lib_dirs, libs

def get_includes():

    p = subprocess.Popen(['pkg-config', '--cflags', 'dip'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout_output, stderr_output = p.communicate()

    if p.returncode != 0 :
        raise Exception('pkg-config terminated with status: %d. stderr' +
                        'follows: %s' % (p.returncode, stderr_output))
    elif stderr_output:
        print stderr_output

    inc_dirs = [flag[2:] for flag in stdout_output.split()
                if flag.startswith('-I')]

    return inc_dirs

def get_defines():

    p = subprocess.Popen(['pkg-config', '--cflags', 'dip'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout_output, stderr_output = p.communicate()

    if p.returncode != 0 :
        raise Exception('pkg-config terminated with status: %d. stderr' +
                        'follows: %s' % (p.returncode, stderr_output))
    elif stderr_output:
        print stderr_output

    defines = [flag[2:] for flag in stdout_output.split()
               if flag.startswith('-D')]

    return defines

def read_file(file_name):

    file_path = join(
        dirname(__file__),
        file_name
        )

    return open(file_path).read()

pkg_config_path = os.environ['PKG_CONFIG_PATH']

if pkg_config_path is None:
    raise Exception('Please set the environment variable PKG_CONFIG_PATH' +
                    'to the location of the COIN pkgconfig directory')

library_dirs, libraries = get_libs()

include_dirs = get_includes()

macros = get_defines()

macros.append(('__DECOMP_LP_CLP__', None))

files = ['DippyDecompAlgo.cpp',
         'DippyDecompApp.cpp',
         'DippySolve.cpp',
         'DippyModule.cpp',
         'DippyPythonUtils.cpp',
         ]

sources = [join('src/dippy', f) for f in files]

modules=[Extension('_dippy', 
                   sources, 
                   libraries=libraries,
                   include_dirs=include_dirs,
                   library_dirs=library_dirs,
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


