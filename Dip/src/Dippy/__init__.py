import sys

if sys.platform == 'linux2':
    try:
        import path
    except ImportError:
        pass

    import ctypes

    ctypes.CDLL('libm.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libblas.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('liblapack.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libz.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libbz2.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libCoinUtils.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libOsi.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libClp.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libOsiClp.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libCgl.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libCbc.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libCbcSolver.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libOsiCbc.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libAlps.so', mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL('libDecomp.so', mode=ctypes.RTLD_GLOBAL)
    
from dippy import *

