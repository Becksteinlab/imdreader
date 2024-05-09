from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Compiler import Options
import numpy

# These are optional
Options.docstrings = True
Options.annotate = False

# Modules to be compiled and include_dirs when necessary
extensions = [
    Extension(
        "*",
        # include_dirs is the list of directories to search for C/C++ header files
        ["*.pyx"], include_dirs=[numpy.get_include()],
    ),
]


# This is the function that is executed
setup(
    name='streamreader', 

    # A list of compiler Directives is available at
    # https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#compiler-directives

    # external to be compiled
    ext_modules = cythonize(extensions, compiler_directives={"language_level": 3, "profile": False}),
)