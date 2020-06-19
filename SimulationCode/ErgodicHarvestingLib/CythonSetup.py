from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

print("Compiling Cython codes for accelerated execution...")
extension = Extension(
    name="cyEntropy",
    sources=["cyEntropy.pyx"],
    extra_compile_args=["-Ofast", "-flto", "-march=native", "-fopenmp"],
    extra_link_args=["-Ofast", "-flto", "-march=native", "-fopenmp"],
    define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
)
setup(
    ext_modules=cythonize(
        extension,
        annotate=False,
        language_level=3,
        nthreads=2,
        quiet=True,
        compiler_directives={
            "boundscheck": False,
            "wraparound": False,
            "initializedcheck": False,
            "cdivision": True,
            "language_level": 3,
        },
    ),
    include_dirs=[numpy.get_include()],
)
print("Success!")
