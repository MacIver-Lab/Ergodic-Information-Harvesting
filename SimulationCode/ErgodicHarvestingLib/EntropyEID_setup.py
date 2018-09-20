from distutils.core import setup
from Cython.Build import cythonize
import numpy
import scipy
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True


setup(
    #ext_modules=cythonize(['EntropyEIDCython.pyx']),
    ext_modules=cythonize(['cyEntropy.pyx', 'cyODE.pyx'], annotate=True),
    include_dirs=[numpy.get_include(), scipy.get_include()]
)
