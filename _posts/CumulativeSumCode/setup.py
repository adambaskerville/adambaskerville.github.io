from distutils.core import setup
from Cython.Build import cythonize

extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp"],
extra_link_args=["-fopenmp"]
setup(
    ext_modules = cythonize("CumulativeSumCythonForLoops.pyx", annotate=True)
)
